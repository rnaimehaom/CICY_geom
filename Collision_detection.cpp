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
    //hout<<"gnp1.vertices[0]="<<gnp1.vertices[0].str()<<endl;
    //hout<<"gnp_new.vertices[0]"<<gnp_new.vertices[0].str()<<endl;
    //hout<<"D0="<<D.str()<<endl;
    
    //Initialize the simplex with the support point of direction D
    Point_3D S = Support_AB(D, gnp1.vertices, gnp_new.vertices);
    simplex.push_back(S);
    
    //Initialize the search direction with the vector from S to the origin
    D = S*(-1);
    //hout<<"D="<<D.str()<<endl;
    
    //If the first vertex is zero, then there is a vertex to vertex touch
    if (Support_AB(D, gnp1.vertices, gnp_new.vertices).length2() < Zero) {
        
        //Touch case
        //Although it could also be penetration when the two GNPs occupy the same space
        //but this will be very unlikely to hapen, so it is ignored
        t_flag = true;
        
        //No penetration but touch
        p_flag =  false;
        
        //Terminate the function
        //hout<<"Case 0: Vertices shared, no GJK needed"<<endl;
        return 1;
    }
    
    while (true) {
        
        //Get a new support point
        Point_3D A = Support_AB(D, gnp1.vertices, gnp_new.vertices);
        //hout<<"A="<<A.str()<<endl;
        
        //Check if there is no intersection
        if (A.dot(D) < Zero) {
            
            //This conditions means there is no intersection because the vector OA
            //goes in the direction opposite to D
            //That is, the origin is beyond A. If the origin is beyond A, then it
            //has to be outside the Minkowski sum as A is the furthest point in direction D
            //thus terminate the function with false
            p_flag = false;
            
            //Terminate the function
            //hout<<"Case GJK 1: Origin is outside the Minkowski sum"<<endl;
            return 1;
        }
        
        //If this part of the code is reached, then update the simplex and check if it
        //contains the origin
        if (Is_origin_in_simplex(simplex, A, D, t_flag)) {
            
            //Add the last vertex to the simplex
            simplex.push_back(A);
            
            //If the simplex contains the origin, set the penetration flag to true
            p_flag = true;
            
            //Terminate the function
            //hout<<"Case GJK 2: Origin is inside the Minkowski sum"<<endl;
            return 1;
        }
        else if (t_flag && Check_actual_touch(gnp1.vertices, gnp_new.vertices, D, simplex)) {
            
            //Set the flag for touching to true as the polyhedrons are touching
            t_flag = true;
            
            //Set the penetration flag to false, as touching is not penetration
            p_flag = false;
            
            //Terminate the function
            //hout<<"Case GJK 3: Touch"<<endl;
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
bool Collision_detection::Is_origin_in_simplex(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, bool &t_flag) {
    
    if (simplex.size() == 1) {
        
        //hout<<"case simplex 1"<<endl;
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
            
            //The origin is closest to vertex A, so the simplex is only A
            simplex[0] = A;
            //hout<<"simplex={A="<<simplex[0].str()<<"}"<<endl;
            
            //The new direction is from A to the origin
            D = AO;
            //hout<<"D=AO="<<D.str()<<endl;
        }
        
        //Return false as either way we cannot say that the origin is inside the simplex
        return false;
    }
    else if (simplex.size() == 2) {
        //Go to the case where the simplex has two elements before (potentially)
        //adding the new vertex A
        return Update_simplex_case2(simplex, A, D, t_flag);
    }
    else if (simplex.size() == 3) {
        //Go to the case where the simplex has three elements before (potentially)
        //adding the new vertex A
        return Update_simplex_case3(simplex, A, D, t_flag);
    }
    else {
        hout<<"There was a problem in the code, simplex has "<<simplex.size()<<" elments"<<endl;
        //Also return true to terminate the function GJK
        return true;
    }
}
//This is the case when the simplex already has two points and the support point A is added
bool Collision_detection::Update_simplex_case2(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, bool &t_flag)
{
    //hout<<"case simplex 2"<<endl;
    
    //Note that simplex = {B, C}
    //Point_3D B = simplex[0];
    //Point_3D C = simplex[1];
    //hout<<"simplex0={B="<<simplex[0].str()<<", C="<<simplex[1].str()<<"}"<<endl;
    
    //Precompute some vectors
    Point_3D AB = simplex[0]-A;
    Point_3D AC = simplex[1]-A;
    Point_3D AO = A*(-1);
    
    //ABC = ACxAB
    Point_3D ABC = AC.cross(AB);
    //hout<<"ABC="<<ABC.str()<<endl;
    
    //Calculate the normal to edge AC, such that it is on
    //the plane of the triangle and goes outside the triangle
    Point_3D N_AC = AC.cross(ABC);
    //hout<<"N_AC="<<N_AC.str()<<endl;
    
    //Check if the origin is outside the triangle, on the side of edge AC
    //or in the plane of the edge
    double dot_N_AC = N_AC.dot(AO);
    //hout<<"N_AC.dot(AO)="<<dot_N_AC<<endl;
    if (dot_N_AC > Zero || abs(dot_N_AC) < Zero) {
        
        //Check if the origin is on the region closest to edge AC
        //hout<<"AC.dot(AO)="<<AC.dot(AO)<<endl;
        if (AC.dot(AO) > 0) {
            
            //Update simplex as {A, C}
            //Just substitute B with A
            simplex[0] = A;
            //hout<<"simplex={A="<<simplex[0].str()<<", C="<<simplex[1].str()<<"}"<<endl;
            
            //Calculate the direction
            D = AC.cross(AO.cross(AC));
            
            //If the direction vector is zero, then the origin is in the edge
            //and probably there is touch
            if (D.length2() < Zero) {
                
                //Set the touch flag to true
                //hout<<"touch 2b 1"<<endl;
                t_flag = true;
                
                //Set the direction as the normal to the edge in the place of the triangle
                D = N_AC;
            }
            
        }
        else {
            
            //The origin is closest to vertex A, so the simplex is only A
            simplex.assign(1, A);
            //hout<<"simplex={A="<<simplex[0].str()<<"}"<<endl;
            
            //The new direction is from A to the origin
            D = AO;
            //hout<<"D=AO="<<D.str()<<endl;
        }
    }
    else {

        //Calculate the normal to edge AB, such that it is on
        //the plane of the triangle and goes outside the triangle
        Point_3D N_AB = ABC.cross(AB);
        //hout<<"N_AB="<<N_AB.str()<<endl;
        
        //Check if the origin is outside the triangle, on the side of edge AB
        //or in the plane of the edge
        double dot_N_AB = N_AB.dot(AO);
        //hout<<"N_AB.dot(AO)="<<dot_N_AB<<endl;
        if (dot_N_AB > Zero || abs(dot_N_AB) < Zero) {
            
            //Check if the origin is on the region closest to edge AB
            //hout<<"AB.dot(AO)="<<AB.dot(AO)<<endl;
            if (AB.dot(AO) > Zero) {
                
                //Update simplex as {B, A}
                //Just substitute C with A
                simplex[1] = A;
                //hout<<"simplex={B="<<simplex[0].str()<<", A="<<simplex[1].str()<<"}"<<endl;
                
                //Calculate the direction
                D = (AB.cross(AO)).cross(AB);
                
                //If the direction vector is zero, then the origin is in the edge
                //and probably there is touch
                if (D.length2() < Zero) {
                    
                    //Set the touch flag to true
                    //hout<<"touch 2b 2"<<endl;
                    t_flag = true;
                    
                    //Set the direction as the normal to the edge in the place of the triangle
                    D = N_AB;
                }
            }
            else {
                
                //The origin is closest to vertex A, so the simplex is only A
                simplex.assign(1, A);
                //hout<<"simplex={A="<<simplex[0].str()<<"}"<<endl;
                
                //The new direction is from A to the origin
                D = AO;
                //hout<<"D=AO="<<D.str()<<endl;
            }
        }
        else {
            
            //The origin is either above or below the triangle, so the simplex
            //will be {B,C,A}
            simplex.push_back(A);
            
            //Set D to have the direction above the triangle
            D = ABC;
            
            //Check if the direction of D needs to be reversed
            double dot_ABC = ABC.dot(AO);
            //hout<<"ABC.dot(AO)="<<dot_ABC<<endl;
            if (abs(dot_ABC) < Zero) {
                //The origin is inside the triangle (on its surface), likely there is touch
                t_flag = true;
            }
            else if (dot_ABC < Zero) {
                //Just reverse the direction of D
                D = D*(-1.0);
            }
        }
    }
    
    //Terminate the function with false, as the origin is not inside the simplex
    //in any of the cases above
    return false;
}
//This is the case when the simplex already has three points and the support point A is added
bool Collision_detection::Update_simplex_case3(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, bool &t_flag)
{
    //hout<<"case simplex 3"<<endl;
    
    //Get a vector from tetrahedron to origin
    Point_3D AO = A*(-1);
    
    //Rename points for simplicity
    //Recall that simplex = {B, C, E}
    Point_3D B = simplex[0];
    Point_3D C = simplex[1];
    Point_3D E = simplex[2];
    //hout<<"simplex0={B="<<simplex[0].str()<<", C="<<simplex[1].str()<<", E="<<simplex[2].str()<<"}"<<endl;
    
    //Get the normals to the faces that go outside the tetrahedron
    Point_3D ACE = Normal_to_face_ouside_tetrahedron(A, C, E, B);
    //hout<<"ACE="<<ACE.str()<<endl;
    Point_3D AEB = Normal_to_face_ouside_tetrahedron(A, E, B, C);
    //hout<<"AEB="<<AEB.str()<<endl;
    Point_3D ABC = Normal_to_face_ouside_tetrahedron(A, B, C, E);
    //hout<<"ABC="<<ABC.str()<<endl;
    
    //Get the dot product ACE.AO
    double dot_ACE = ACE.dot(AO);
    //hout<<"ACE.dot(AO)="<<dot_ACE<<endl;
    
    //Check if the origin is in the face ACE
    if (abs(dot_ACE) < Zero) {
        
        //The origin is in the face ACE, probably there is touch
        //hout<<"touch in case 3 ACE"<<endl;
        t_flag = true;
        
        //Set the simplex as the face ACE
        simplex[0] = A;
        
        //Terminate with false
        return false;
    }
    
    //Check if the origin is above the plane of face ACE
    if (dot_ACE > Zero) {
        
        //Since origin is above the plane of face ACE, go to case 2 for tetrahedron
        //with simplex = {C,E}, A and B
        simplex[0] = C;
        simplex[1] = E;
        simplex.pop_back();
        //hout<<"simplex to 2 for tetrahedron={C="<<simplex[0].str()<<", E="<<simplex[1].str()<<"}"<<endl;
        
        return Update_simplex_case2_tetrahedron(simplex, A, D, B, ACE, ABC, AEB, t_flag);
    }
    else {
        
        //Get the dot product ABC.AO
        double dot_AEB = AEB.dot(AO);
        //hout<<"AEB.dot(AO)="<<dot_AEB<<endl;
        
        //Check if the origin is in the face AEB
        if (abs(dot_AEB) < Zero) {
            
            //The origin is in the face AEB, probably there is touch
            //hout<<"touch in case 3 AEB"<<endl;
            t_flag = true;
            
            //Set the simplex as the face AEB
            simplex[1] = A;
            
            //Terminate with false
            return false;
        }
        
        //Check if the origin is above the plane of face AEB
        if (dot_AEB > Zero) {
            
            //Since origin is above the plane of face AEB, go to case 2 for tetrahedron
            //with simplex = {E, B}, A and C
            simplex[0] = E;
            simplex[1] = B;
            simplex.pop_back();
            //hout<<"simplex to 2 for tetrahedron={E="<<simplex[0].str()<<", B="<<simplex[1].str()<<"}"<<endl;
            
            return Update_simplex_case2_tetrahedron(simplex, A, D, C, AEB, ACE, ABC, t_flag);
        }
        else {
            
            //Get the dot product ABC.AO
            double dot_ABC = ABC.dot(AO);
            //hout<<"ABC.dot(AO)="<<dot_ABC<<endl;
            
            //Check if the origin is in the face ABC
            if (abs(dot_ABC) < Zero) {
                
                //The origin is in the face ABC, probably there is touch
                //hout<<"touch in case 3 ABC"<<endl;
                t_flag = true;
                
                //Set the simplex as the face ABC
                simplex[2] = A;
                
                //Terminate with false
                return false;
            }
            
            //Check if the origin is above the plane of face ABC
            if (dot_ABC > Zero) {
                
                //Since origin is above the plane of face ABC, go to case 2
                //with simplex = {B, C}, A and E
                simplex.pop_back();
                //hout<<"simplex to 2 for tetrahedron={B="<<simplex[0].str()<<", C="<<simplex[1].str()<<"}"<<endl;
                
                return Update_simplex_case2_tetrahedron(simplex, A, D, E, ABC, AEB, ACE, t_flag);
            }
            else {
                
                //The origin is inside the Minkowski sum, so there is nothing to do
                //hout<<"Origin is inside the simplex (Case 3)"<<endl;
                
                //Return true as the origin is inside the simplex
                return true;
            }
        }
    }
}

bool Collision_detection::Update_simplex_case2_tetrahedron(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, Point_3D &E, Point_3D &ABC, Point_3D &AEB, Point_3D &ACE, bool &t_flag)
{
    //hout<<"case simplex 2 tetrahedron"<<endl;
    
    //Note that simplex = {B, C}
    //hout<<"simplex0={B="<<simplex[0].str()<<", C="<<simplex[1].str()<<"}"<<endl;
    
    //Precompute some vectors
    //Note that simplex = {B, C}
    Point_3D AB = simplex[0]-A;
    Point_3D AC = simplex[1]-A;
    Point_3D AO = A*(-1);
    
    //ABC = ACxAB
    //hout<<"ABC="<<ABC.str()<<endl;
    
    //Calculate the normal to edge AC, such that it is on
    //the plane of the triangle and goes outside the triangle
    Point_3D N_AC = AC.cross(ABC);
    //Make sure it goes outside te triangle
    if (N_AC.dot(AB) > Zero) {
        N_AC = N_AC*(-1);
    }
    //hout<<"N_AC="<<N_AC.str()<<endl;
    
    //Check if the origin is outside the triangle, on the side of edge AC
    //or in the plane of the edge
    double dot_N_AC = N_AC.dot(AO);
    //hout<<"N_AC.dot(AO)="<<dot_N_AC<<endl;
    if (dot_N_AC > Zero || abs(dot_N_AC) < Zero) {
        
        //Check if the origin is outside the triangle, on the side of edge AC
        //hout<<"AC.dot(AO)="<<AC.dot(AO)<<endl;
        if (AC.dot(AO) > 0) {
            
            //Get the Normal to edge AC on the plane of face ACE
            Point_3D N_AC_ACE = AC.cross(ACE);
            //Make sure it goes outside the triangle ACE
            if (N_AC_ACE.dot(E-A) > Zero) {
                N_AC_ACE = N_AC_ACE*(-1);
            }
            //hout<<"N_AC_ACE="<<N_AC_ACE.str()<<endl;
            
            //Check if origin is closest to edge AC or face shared with edge AC (face ACE)
            //hout<<"N_AC_ACE.dot(AO)="<<N_AC_ACE.dot(AO)<<endl;
            if (N_AC_ACE.dot(AO) > Zero) {
                
                //Update simplex as {A, C}
                //Just substitute B with A
                simplex[0] = A;
                //hout<<"simplex={A="<<simplex[0].str()<<", C="<<simplex[1].str()<<"}"<<endl;
                
                //Calculate the direction
                D = AC.cross(AO.cross(AC));
                
                //If the direction vector is zero, then the origin is in the edge
                //and probably there is touch
                if (D.length2() < Zero) {
                    
                    //Set the touch flag to true
                    //hout<<"touch 2b 1"<<endl;
                    t_flag = true;
                    
                    //Set the direction as the normal to the edge in the place of the triangle
                    D = N_AC;
                }
            }
            else {
                //The origin is actually closer to face ACE
                
                //Update simplex as {A, C, E}
                simplex[0] = A;
                simplex.push_back(E);
                //hout<<"simplex={A="<<simplex[0].str()<<", C="<<simplex[1].str()<<", E="<<simplex[2].str()<<"}"<<endl;
                
                //The direction is the normal to face ACE
                D = ACE;
            }
        }
        else {
            
            //The origin is closest to vertex A, so the simplex is only A
            simplex.assign(1, A);
            //hout<<"simplex={A="<<simplex[0].str()<<"}"<<endl;
            
            //The new direction is from A to the origin
            D = AO;
            //hout<<"D=AO="<<D.str()<<endl;
        }
    }
    else {

        //Calculate the normal to edge AB, such that it is on
        //the plane of the triangle and goes outside the triangle
        Point_3D N_AB = ABC.cross(AB);
        //Make sure it goes outside te triangle
        if (N_AB.dot(AC) > Zero) {
            N_AB = N_AB*(-1);
        }
        //hout<<"N_AB="<<N_AB.str()<<endl;
        
        //Check if the origin is outside the triangle, on the side of edge AB
        //or in the plane of the edge
        double dot_N_AB = N_AB.dot(AO);
        //hout<<"N_AB.dot(AO)="<<dot_N_AB<<endl;
        if (dot_N_AB > Zero || abs(dot_N_AB) < Zero) {
            
            //Check if the origin is on the region closest to edge AB
            //hout<<"AB.dot(AO)="<<AB.dot(AO)<<endl;
            if (AB.dot(AO) > Zero) {
                
                //Get the Normal to edge AB on the plane of face AEB
                Point_3D N_AB_AEB = AEB.cross(AB);
                //hout<<"N_AB_AEB="<<N_AB_AEB.str()<<endl;
                //Make sure it goes outside the triangle AEB
                if (N_AB_AEB.dot(E-A) > Zero) {
                    N_AB_AEB = N_AB_AEB*(-1);
                }
                
                //Check if origin is closest to edge AB or face shared with edge AB (face AEB)
                //hout<<"N_AB_AEB.dot(AO)="<<N_AB_AEB.dot(AO)<<endl;
                if (N_AB_AEB.dot(AO) > Zero) {
                    //Update simplex as {B, A}
                    //Just substitute C with A
                    simplex[1] = A;
                    //hout<<"simplex={B="<<simplex[0].str()<<", A="<<simplex[1].str()<<"}"<<endl;
                    
                    //Calculate the direction
                    D = (AB.cross(AO)).cross(AB);
                    
                    //If the direction vector is zero, then the origin is in the edge
                    //and probably there is touch
                    if (D.length2() < Zero) {
                        
                        //Set the touch flag to true
                        //hout<<"touch 2b 2"<<endl;
                        t_flag = true;
                        
                        //Set the direction as the normal to the edge in the place of the triangle
                        D = N_AB;
                    }
                }
                else {
                    //The origin is actually closer to face AEB
                    
                    //Update simplex as {B, A, E}
                    simplex[1] = A;
                    simplex.push_back(E);
                    //hout<<"simplex={B="<<simplex[0].str()<<", A="<<simplex[1].str()<<", E="<<simplex[2].str()<<"}"<<endl;
                    
                    //The direction is the normal to face AEB
                    D = AEB;
                }
            }
            else {
                
                //The origin is closest to vertex A, so the simplex is only A
                simplex.assign(1, A);
                //hout<<"simplex={A="<<simplex[0].str()<<"}"<<endl;
                
                //The new direction is from A to the origin
                D = AO;
                //hout<<"D=AO="<<D.str()<<endl;
            }
        }
        else {
            
            //The origin is either above or below the triangle, so the simplex
            //will be {B,C,A}
            simplex.push_back(A);
            //hout<<"simplex={B="<<simplex[0].str()<<", C="<<simplex[1].str()<<", A="<<simplex[2].str()<<"}"<<endl;
            
            //Set D to have the direction above the triangle
            D = ABC;
            
            //Check if the direction of D needs to be reversed
            double dot_ABC = ABC.dot(AO);
            //hout<<"ABC.dot(AO)="<<dot_ABC<<endl;
            if (abs(dot_ABC) < Zero) {
                //The origin is inside the triangle (on its surface), likely there is touch
                //hout<<"Touch in case 2"<<endl;
                t_flag = true;
            }
            else if (dot_ABC < Zero) {
                //Just reverse the direction of D
                //hout<<"Reverse D"<<endl;
                D = D*(-1.0);
            }
        }
    }
    
    //Terminate the function with false, as the origin is not inside the simplex
    //in any of the cases above
    return false;
}
//When the t_flag is set to true, the origin is in either an edge or a face of the simplex.
//This function helps finding if the origin is inside the Minkowski sum or on its edge or face.
//If the origin is in the edge or face of the Minkowski sum, then there is touch. Othwewise, the
//origin is in an edge or face of the simplex, but this edge or face is inside the Minkowski sum.
bool Collision_detection::Check_actual_touch(const Point_3D verticesA[], const Point_3D verticesB[], const Point_3D &D, vector<Point_3D> &simplex)
{
    
    //Get a support point in the direction of D
    Point_3D S = Support_AB(D, verticesA, verticesB);
    
    if (simplex.size() == 3) {
        
        //Check if S is in the plane of the simplex
        //This happens if a vector from a vertex of the simplex to S has zero dot product
        //With the normal to the plane
        Point_3D N = (simplex[1] - simplex[0]).cross(simplex[2] - simplex[0]);
        if ( abs((S - simplex[0]).dot(N)) < Zero) {
            
            //There is actual touch as the simplex is a face of the Minkowski sum
            return true;
        }
        
    }
    else if (simplex.size() == 2) {
        
        //Check if S is still in the edge defined by the simplex
        //This happens if the vector from a vertex of the edge to S has zero dot with D
        if (abs( (S - simplex[0]).dot(D) )) {
            
            //There is actual touch as the simplex is an edge of the Minkowski sum
            return true;
        }
    }
    return false;
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
    //Face BCD is ignored because by construction, it cannot be the closest face to the origin
    //as it does not contain A. By construction A is the closest point to the origin.
    //So a face that does not contain A cannot be closest to the origin
    
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
    
    //Calculate the dot products of the normals to the faces with the vector AO, and check the signs
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
                hout<<"Error in Find_simplex_closest_to_origin (1). Vector idxs has an invalid combination:"<<endl;
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
            hout<<"Error in Find_simplex_closest_to_origin (2). Vector idxs has an invalid combination:"<<endl;
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
