//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
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
    D = S*(-1.0);
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
        //hout<<"A="<<A.str()<<endl<<"D="<<D.str()<<endl;
        
        //Check if there is no penetration, i.e., the origin is furthest than th support point
        if (A.dot(D) < Zero) {
            
            //This conditions means there is no penetration because the vector OA
            //goes in the direction opposite to D
            //That is, the origin is beyond A. If the origin is beyond A, then it
            //has to be outside the Minkowski sum as A is the furthest point in direction D
            //thus set the penetration flag to false
            p_flag = false;
            
            //This last update is not actually needed since the update of the simplex
            //using the function ignores certain simplices of this simplex. In this final
            //Step all simplices of the simplex need to be check when searching for
            //the final simplex that is closest to the origin. Because of this, the simplex
            //that is returned when making this update is not always the closest to the
            //origin. Thus, this update is not done anymore. For now it is clear to me that
            //no update is needed and the simplex found by the GJK is the final simplex and
            //support point A is not part of that simplex (unless it already is part of the
            //final simplex)
            //However, I leave the commented code in case some variation of this idea is needed
            /* /Make a last update to the simplex if needed
            if(Is_point_in_simplex(simplex, A)) {
                
                //hout<<"A already in simplex"<<endl;
                
                //Then take the last point in simplex to be A
                A = simplex.back();
                
                //Remove the last point in the simplex so that
                //one last update to the simplex can be done
                simplex.pop_back();
            }
            //else {hout<<"A not in simplex"<<endl;}
            
            //Get the closest simplex to the origin (one last time)
            //There is no need to store the output of the function is_origin_in_simplex
            //hout<<"Last update of the simplex"<<endl;
            Is_origin_in_simplex(simplex, A, D, t_flag);*/
            
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
    return (Support_map(D, verticesA) - Support_map(D*(-1.0), verticesB) );
}
//Function for the suport point in a given direction
Point_3D Collision_detection::Support_map(const Point_3D &D, const Point_3D vertices[])
{
    //Initialize distance with first index
    double dist = D.dot(vertices[0]);
    
    //Initialize index of maximum dot product with 0 (first index)
    int idx = 0;
    //hout.precision(17);
    //hout << "i=" << idx << " D=" << D.str() << " v[i]=" << vertices[0].str() << " dist=" << dist << endl;
    
    //Find the maximum distance in the rest of the vertices
    //The array of vertices has a fixed size of 8, so that value is used here
    //It will need to change if other types of polyhedrons are used
    for (int i = 1; i < 8; i++) {
        
        //Calculate the new distance
        double new_dist = D.dot(vertices[i]);
        //hout << "i=" << i << " D=" << D.str() << " v[i]=" << vertices[i].str() << " dist=" << new_dist <<" new_dist-dist="<< new_dist - dist << endl;
        
        //Check if update is needed
        if (new_dist - dist > Zero) {
            //Update variables
            dist = new_dist;
            idx = i;
        }
    }
    //hout << "support idx=" << idx << " dist="<< D.dot(vertices[idx]) << endl;

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
        Point_3D AO = A*(-1.0);
        
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
    Point_3D AO = A*(-1.0);
    
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
        if (AC.dot(AO) > Zero) {
            
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
    Point_3D AO = A*(-1.0);
    
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
//This function finds the normal to a tetrahedron face that goes outside the tetrahedron
Point_3D Collision_detection::Normal_to_face_ouside_tetrahedron(const Point_3D &A, const Point_3D &B, const Point_3D &C, const Point_3D &not_in_face)
{
    //Get a normal to face ABC
    Point_3D ABC = (B - A).cross(C - A);
    
    //Make sure ABC goes outside the tetrahedron, e.g., in the oposite direction of AD,
    //where D is not_in_face
    Point_3D AD = not_in_face - A;
    
    //Calculate dot product
    double dot_tmp = ABC.dot(AD);
    if (dot_tmp > Zero) {
        //Reverse the direction is ABC and AD go in the same direction, i.e.,
        //ABC goes inside the tetrahedron
        ABC = ABC*(-1.0);
    }
    
    return ABC;
}
//Find the simplex closest to the origin in the case that A and simplex are a face of a tetrahedron
bool Collision_detection::Update_simplex_case2_tetrahedron(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, Point_3D &E, Point_3D &ABC, Point_3D &AEB, Point_3D &ACE, bool &t_flag)
{
    //hout<<"case simplex 2 tetrahedron"<<endl;
    
    //Note that simplex = {B, C}
    //hout<<"simplex0={B="<<simplex[0].str()<<", C="<<simplex[1].str()<<"}"<<endl;
    
    //Precompute some vectors
    //Note that simplex = {B, C}
    Point_3D AB = simplex[0]-A;
    Point_3D AC = simplex[1]-A;
    Point_3D AO = A*(-1.0);
    
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
        if (AC.dot(AO) > Zero) {
            
            //Get the Normal to edge AC on the plane of face ACE
            Point_3D N_AC_ACE = AC.cross(ACE);
            //Make sure it goes outside the triangle ACE
            if (N_AC_ACE.dot(E-A) > Zero) {
                N_AC_ACE = N_AC_ACE*(-1.0);
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
            N_AB = N_AB*(-1.0);
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
                    N_AB_AEB = N_AB_AEB*(-1.0);
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
//This function checks is the support point A is already in the simplex
bool Collision_detection::Is_point_in_simplex(const vector<Point_3D> &simplex, const Point_3D &A)
{
    
    for (size_t i = 0; i < simplex.size(); i++) {
        //hout << "simplex[" << i << "]=" << simplex[i].str() << endl;
        if (simplex[i].squared_distance_to(A) < Zero) {
            //hout<<"A is already in the simplex"<<endl;
            return true;
        }
    }
    return false;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//GJK for distance calculation
int Collision_detection::GJK_distance(const GNP& gnp1, const GNP& gnp2, vector<Point_3D>& simplex, double& dist, Point_3D& N, bool& p_flag)
{
    //hout << "=========================================" << endl;
    //hout << "GNP1=" << gnp1.flag << " GNP2=" << gnp2.flag << endl;
    //for (size_t i = 0; i < simplex.size(); i++)
    //    hout << "simplex[" << i << "]=" << simplex[i].str() << endl;

    //Check if the simplex is empty
    if (simplex.empty()) {

        //The simplex is empty, this is the case when no previous boolean GJK has been used
        //Initialize the simplex with a vertex form the Minkowski sum
        //hout << "simplex is empty" << endl;
        simplex.push_back(gnp1.vertices[0] - gnp2.vertices[0]);
    }

    //Calculate initial distance from simplex to origin
    if (!Distance_from_simplex_to_origin(simplex, dist)) {
        hout << "Error in GJK_distance when calling Distance_from_simplex_to_origin (0)" << endl;
        return 0;
    }

    //Iteration count
    int it = 1;

    //Iterate up to the maximum number of iterations
    while (it <= 50) {

        //hout << "===it " << it << endl; 

        //Calculate the new direction
        Point_3D D;
        if (!Direction_from_simplex_to_origin(simplex, D)) {
            hout << "Error in GJK_distance when callin Direction_from_simplex_to_origin" << endl;
            return 0;
        }

        //Get a support point in the Minkowski sum
        Point_3D A = Support_AB(D, gnp1.vertices, gnp2.vertices);
        //hout << "A=" << A.str() << " D=" << D.str() << " A.dot(D)=" << A.dot(D) << endl;

        //Check if A is already in the simplex
        //hout << "is_point_in_simplex" << endl;
        int s_flag = Is_point_contained_in_simplex(simplex, A);
        if (s_flag == 1) {

            //hout << "dist=" << dist << " |A.dot(D)|=" << abs(A.dot(D)) << endl;
            //Check if simplex and support point are at the same distance to the origin
            //if (abs(dist - abs(A.dot(D))) < Zero) {

                //hout << "Case 2: Simplex cannot be updated as support point is contained in simplex" << endl;
                //Distance from simplex to origin is not calculated again because it has
                //not changed as the simplex has not changed since this distance was calculated
                //Get the direction
                N = D;


                return 1;
            //}
            //If the simplex and support point are not at the same distance, 
            //then continue with the next iteration
            //hout << "Support point in simplex but termination criteria not met" << endl;

        }
        else if (s_flag == 0) {

            //Add support point A to the simplex if not contained in it
            simplex.push_back(A);
            //hout << "updated simplex.size=" << simplex.size() << endl;
        }
        else if (s_flag == -1) {
            hout << "Error in GJK_distance_from_simplex_revised when calling is_point_contained_in_simplex" << endl;
            return 0;
        }

        //Find the simplex within the vector simplex that is closest to the origin
        int vr = Find_voroni_region(simplex);

        //Check if there is penetration
        if (vr == 0) {
            //Set the penetration flag to true and terminate the function
            p_flag = true;
            return 1;
        }
        //Check for error
        if (vr == -1) {
            hout << "Error in GJK_distance when calling Find_voroni_region" << endl;
            return 0;
        }
        //hout << "voroni simplex.size=" << simplex.size() << endl;

        //Calculate the distance from the updated simplex to the origin
        if (!Distance_from_simplex_to_origin(simplex, dist)) {
            hout << "Error in GJK_distance when calling Distance_from_simplex_to_origin (1)" << endl;
            return 0;
        }
        //hout << "dist=" << dist << endl << endl;

        /*if (dist < Zero) {
            hout << "Touch" << endl;
        }*/

        //Check if support point is closer to origin than simplex
        if (abs(dist - abs(A.dot(D))) < Zero) {

            //hout << "Case 1: Simplex is as close to origin as support point" << endl;
            //hout << "dist=" << dist << " abs(A.dot(D))=" << abs(A.dot(D)) << " A.length=" << A.length() << " A=" << A.str() << endl;
            //Calculate direction
            N = D;

            return 1;
        }

        //Increase the number of iterations
        it++;
    }

    hout << "Error in GJK_distance: Maximum number of iterations reached." << endl;
    return 0;
}
//This function calculates the directino vector from a simplex towards the origin
int Collision_detection::Direction_from_simplex_to_origin(const vector<Point_3D>& simplex, Point_3D &N)
{
    //Check the size of the simplex
    if (simplex.size() == 1) {

        //The direction vector is AO
        N = (simplex[0] * (-1.0)).unit();
    }
    else if (simplex.size() == 2) {

        //Calculate some vectors
        Point_3D AO = simplex[0] * (-1.0);
        //AB = B - A
        Point_3D AB = simplex[1] - simplex[0];

        //The direction vector is the normal from AB going towards O
        //N = (ABxAO)xAB
        N = ((AB.cross(AO)).cross(AB)).unit();
    }
    //simplex.size() == 3
    else if (simplex.size() == 3) {

        //Get a unit normal to the plane defined by the three points
        //This normal is the direction vector
        N = (simplex[1] - simplex[0]).cross(simplex[2] - simplex[0]);
        N.make_unit();

        //Check that it goes to the origin
        if (N.dot(simplex[0]) > Zero) {
            N = N * (-1.0);
        }
    }
    else {
        hout << "Error in Direction_from_simplex_to_origin. Invalid size of simplex to calculate direction from simplex towards origin. Elements in simplex: " << simplex.size() << endl;
        return 0;
    }

    return 1;
}
//This function determines if a point is conteined in the simplex
//This is done by calculating the distance from the point to the simplex, in this way not
//only it is considered the case when one vertex in the simples is the point A, but also
//the case when the point A is in the same plane or line as the simplex
int Collision_detection::Is_point_contained_in_simplex(const vector<Point_3D>& simplex, const Point_3D& A)
{
    if (simplex.size() == 1) {

        //Simplex is a single point, so calculate the distance between points
        if (simplex[0].distance_to(A) < Zero) {

            //hout << "Dist from simplex1=" << simplex[0].distance_to(A) << endl;
            //Points are practically the same, so terminate with true
            return 1;
        }
    }
    else if (simplex.size() == 2) {

        //Get the normal to the plane defined by the two vertices in the simplex and A
        //S0S1xS0A
        Point_3D S0A = A - simplex[0];
        Point_3D S0S1 = simplex[1] - simplex[0];
        Point_3D N_SA = S0S1.cross(S0A);

        //Get the normal to the simplex that goes towards A
        Point_3D N = N_SA.cross(S0S1);

        //Distance from simplex to point A is the dot product of S0A with N (where N is a unit vector)
        //hout << "Dist from simplex2=" << N.dot(S0A) << endl;
        if (N.dot(S0A) < Zero) {

            //Distance from simplex to point A is zero, so A is contained in the simplex
            return 1;
        }
    }
    else if (simplex.size() == 3) {

        //Get the normal to the plane defined by the simplex
        Point_3D N = (simplex[1] - simplex[0]).cross(simplex[2] - simplex[0]);

        //Distance from simplex to point A is the dot product of S0A with N
        //hout << "Dist from simplex3=" << abs(N.dot(A - simplex[0])) << endl;
        if (abs(N.dot(A - simplex[0])) < Zero) {

            //Distance from simplex to point A is zero, so A is contained in the simplex
            return 1;
        }
    }
    else {
        //hout << "Error in is_point_contained_in_simplex. Invalid size of simplex. Size: " << simplex.size() << endl;
        return -1;
    }

    //If this part of the code is reached, then point A is not contained in the simplex
    return 0;
}
//This function calculates the distance from 
int Collision_detection::Distance_from_simplex_to_origin(vector<Point_3D>& simplex, double &dist) 
{
    //Check the size of the simplex
    if (simplex.size() == 1) {

        //hout << "simplex.size=1" << endl;
        //Just calculate the length of the vector OA
        dist = simplex[0].length();
    }
    else if (simplex.size() == 2) {

        //hout << "simplex.size=2" << endl;
        //Calculate some vectors
        Point_3D AO = simplex[0] * (-1.0);
        //AB = B - A
        Point_3D AB = simplex[1] - simplex[0];

        //dp^2 = (AO.dot(AB))^2/||AB||^2
        double dot = AO.dot(AB);
        double dp_2 = dot * dot / AB.length2();

        dist = sqrt(AO.length2() - dp_2);
    }
    //simplex.size() == 3
    else if (simplex.size() == 3) {

        //hout << "simplex.size=3" << endl;
        //Get a unit normal to the plane defined by the three points
        //This normal is the direction vector
        Point_3D N = (simplex[1] - simplex[0]).cross(simplex[2] - simplex[0]);
        N.make_unit();

        dist = abs(N.dot(simplex[0]));
    }
    else {
        hout<<"Error in Distance_from_simplex_to_origin. Invalid size of simplex to calculate distance from simplex to origin. Elements in simplex: " << simplex.size() << endl;
        return 0;
    }

    return 1;
}
//Find the sub-simplex within the simplex that is closest to the origin
//This is equivalent to find the Voroni region of the simplex in which the origin lies
int Collision_detection::Find_voroni_region(vector<Point_3D>& simplex) 
{
    if (simplex.size() == 2) {
        //hout << "Voroni case simplex.size=2" << endl;

        //Calculate a vector along the simplex
        Point_3D S0S1 = simplex[1] - simplex[0];

        //Calculate some dot products
        double dot_s0 = -simplex[0].dot(S0S1);
        double dot_s1 = simplex[1].dot(S0S1);

        //Find the Voroni region
        if (dot_s0 > Zero && dot_s1 > Zero) {
            //hout << "Edge itself is closest to the origin" << endl;
            //Simplex is the closest to the origin, so terminate the function
            return 1;
        }
        else if (dot_s1 < Zero) {
            //Given that the the edge is not closest to the origin,
            //if dot product is negative then vertex simplex[1]
            //is closest to the origin
            //hout << "Vertex 1 from edge" << endl;
            simplex[0] = simplex[1];

            //After leaving the if-statement the second vertex will be removed
        }

        //The remaining option is that vertex simplex[0] is closest
        //to the origin, which case only remove the second vertex
        simplex.pop_back();
        //hout << "Vertex 0 from edge" << endl;
    }
    else if (simplex.size() == 3) {

        //hout << "Voroni case simplex.size=3" << endl;
        return Voroni_region_case3(simplex);
    }
    else if (simplex.size() == 4) {

        //hout << "Voroni case simplex.size=4" << endl;
        return Voroni_region_case4(simplex);
    }
    else if (simplex.size() != 1) {
        hout << "Error. Simplex size is not in the interval [1, 4]. Size: " << simplex.size() << endl;
        return -1;
    }

    return 1;
}
//This fucntion finds the sub-simplex within the simplex that is closest to the origin when 
//the simples has size 3
//This is equivalent to find the Voroni region of the simplex in which the origin lies for
//the case when the simples has size 3
int Collision_detection::Voroni_region_case3(vector<Point_3D>& simplex) 
{
    //hout << "simplex0={S0=" << simplex[0].str() << ", S1=" << simplex[1].str() << ", S2=" << simplex[2].str() << "}" << endl;

    //=============================
    //Vertices
    //hout << endl << "Vertices - -" << endl;

    //Precompute some vectors
    Point_3D S0S1 = simplex[1] - simplex[0];
    Point_3D S0S2 = simplex[2] - simplex[0];

    double s0_s0s1 = -simplex[0].dot(S0S1);
    double s0_s0s2 = -simplex[0].dot(S0S2);

    //Is origin in Voroni region of veretex 0
    if (s0_s0s2 < Zero && s0_s0s1 < Zero) {
        //Update simplex to contain only vertex 0
        //hout << "Vertex 0 from triangle" << endl;
        simplex.pop_back();
        simplex.pop_back();

        //Terminate the function
        return 1;
    }

    //Precompute a vector
    Point_3D S1S2 = simplex[2] - simplex[1];

    //Precompute some dot products
    //This dot product is equivalent to -S1.dot(S1S0):
    double s1_s1s0 = simplex[1].dot(S0S1);
    double s1_s1s2 = -simplex[1].dot(S1S2);

    //Is origin in Voroni region of veretex 1
    if (s1_s1s0 < Zero && s1_s1s2 < Zero) {
        //Update simplex to contain only vertex 1
        //hout << "Vertex 1 from triangle" << endl;
        simplex[0] = simplex[1];
        simplex.pop_back();
        simplex.pop_back();

        //Terminate the function
        return 1;
    }

    //Precompute some dot products
    //This dot product is equivalent to -S2.dot(S2S0):
    double s2_s2s0 = simplex[2].dot(S0S2);
    //This dot product is equivalent to -S2.dot(S2S1):
    double s2_s2s1 = simplex[2].dot(S1S2);

    //Is origin in Voroni region of veretex 2
    if (s2_s2s0 < Zero && s2_s2s1 < Zero) {
        //Update simplex to contain only vertex 2
        //hout << "Vertex 2 from triangle" << endl;
        simplex[0] = simplex[2];
        simplex.pop_back();
        simplex.pop_back();

        //Terminate the function
        return 1;
    }

    //=============================
    //Edges
    //hout << endl << "Edges + + +" << endl;

    //Vector normal to the triangle
    Point_3D NT = S0S1.cross(S0S2);

    //Compute normal to edge S0S1
    Point_3D N01 = S0S1.cross(NT);
    //Make sure normal goes outside the triangle
    if (N01.dot(S1S2) > Zero) {
        N01 = N01 * (-1.0);
    }

    //Is origin in Voroni region of edge S0S1
    if (-simplex[0].dot(N01) > Zero && s0_s0s1 > Zero && s1_s1s0 > Zero) {
        //Update simplex to contain only edge S0S1
        //hout << "Edge 01 from triangle" << endl;
        simplex.pop_back();

        //Terminate the function
        return 1;
    }

    //Compute normal to edge S0S2
    Point_3D N02 = S0S2.cross(NT);
    //Make sure normal goes outside the triangle
    if (N02.dot(S0S1) > Zero) {
        N02 = N02 * (-1.0);
    }

    //Is origin in Voroni region of edge S0S2
    if (-simplex[0].dot(N02) > Zero && s0_s0s2 > Zero && s2_s2s0 > Zero) {
        //Update simplex to contain only edge S0S2
        //hout << "Edge 02 from triangle" << endl;
        simplex[1] = simplex[2];
        simplex.pop_back();

        //Terminate the function
        return 1;
    }

    //Compute normal to edge S1S2
    Point_3D N12 = S1S2.cross(NT);
    //Make sure normal goes outside the triangle
    if (N12.dot(S0S1) < Zero) { //This is equivalent to if (N12.dot(S1S0) > Zero) {
        N12 = N12 * (-1.0);
    }

    //Is origin in Voroni region of edge S1S2
    if (-simplex[1].dot(N12) > Zero && s1_s1s2 > Zero && s2_s2s1 > Zero) {
        //Update simplex to contain only edge S1S2
        //hout << "Edge 12 from triangle" << endl;
        simplex[0] = simplex[1];
        simplex[1] = simplex[2];
        simplex.pop_back();

        //Terminate the function
        return 1;
    }

    //If the function reaches this point, then the simplex closest to the origin
    //is the original simplex

    return 1;
}
//This fucntion finds the sub-simplex within the simplex that is closest to the origin when 
//the simples has size 4
//This is equivalent to find the Voroni region of the simplex in which the origin lies for
//the case when the simples has size 4
int Collision_detection::Voroni_region_case4(vector<Point_3D>& simplex) 
{

    //hout << "case simplex.size=4" << endl;
    //hout << "simplex0={S0=" << simplex[0].str() << ", S1=" << simplex[1].str() << ", S2=" << simplex[2].str() << ", S3=" << simplex[3].str() << "}" << endl;

    //Precompute some vectors
    Point_3D S0S1 = simplex[1] - simplex[0];
    Point_3D S0S2 = simplex[2] - simplex[0];
    Point_3D S0S3 = simplex[3] - simplex[0];

    //Precompute some dot products
    double s0_s0s1 = -simplex[0].dot(S0S1);
    double s0_s0s2 = -simplex[0].dot(S0S2);
    double s0_s0s3 = -simplex[0].dot(S0S3);

    //=============================
    //Vertices
    //hout << "Vertices - - -" << endl;

    //=============================
    //Is origin in Voroni region of veretex 0
    //hout << "Vertex 0: s0_s0s2=" << s0_s0s2 << " s0_s0s1=" << s0_s0s1 << " s0_s0s3=" << s0_s0s3 << endl;
    if (s0_s0s2 < Zero && s0_s0s1 < Zero && s0_s0s3 < Zero) {

        //Update simplex to cotain only vertex 0
        //hout << "Vertex 0 from tetrahedron" << endl;
        if (!Reduce_simplex_4_to_1(0, simplex)) {
            hout << "Error when reducing simplex from 4 to 1. Vertex 0." << endl;
            return -1;
        }

        //Terminate the function
        return 1;
    }

    //Precompute some vectors
    Point_3D S1S2 = simplex[2] - simplex[1];
    Point_3D S1S3 = simplex[3] - simplex[1];

    //Precompute some dot products
    //This dot product is equivalent to -S1.dot(S1S0):
    double s1_s1s0 = simplex[1].dot(S0S1);
    double s1_s1s2 = -simplex[1].dot(S1S2);
    double s1_s1s3 = -simplex[1].dot(S1S3);

    //=============================
    //Is origin in Voroni region of veretex 1
    //hout << "Vertex 1: s1_s1s0=" << s1_s1s0 << " s1_s1s2=" << s1_s1s2 << " s1_s1s3=" << s1_s1s3 << endl;
    if (s1_s1s0 < Zero && s1_s1s2 < Zero && s1_s1s3 < Zero) {

        //Update simplex to cotain only vertex 1
        //hout << "Vertex 1 from tetrahedron" << endl;
        if (!Reduce_simplex_4_to_1(1, simplex)) {
            hout << "Error when reducing simplex from 4 to 1. Vertex 1." << endl;
            return -1;
        }

        //Terminate the function
        return 1;
    }

    //Precompute a vector
    Point_3D S2S3 = simplex[3] - simplex[2];

    //Precompute some dot products
    //This dot product is equivalent to -S2.dot(S2S0):
    double s2_s2s0 = simplex[2].dot(S0S2);
    //This dot product is equivalent to -S2.dot(S2S1):
    double s2_s2s1 = simplex[2].dot(S1S2);
    double s2_s2s3 = -simplex[2].dot(S2S3);

    //=============================
    //Is origin in Voroni region of veretex 2
    //hout << "Vertex 2: s2_s2s0=" << s2_s2s0 << " s2_s2s1=" << s2_s2s1 << " s2_s2s3=" << s2_s2s3 << endl;
    if (s2_s2s0 < Zero && s2_s2s1 < Zero && s2_s2s3 < Zero) {

        //Update simplex to cotain only vertex 2
        //hout << "Vertex 2 from tetrahedron" << endl;
        if (!Reduce_simplex_4_to_1(2, simplex)) {
            hout << "Error when reducing simplex from 4 to 1. Vertex 2." << endl;
            return -1;
        }

        //Terminate the function
        return 1;
    }

    //Precompute some dot products
    //This dot product is equivalent to -S3.dot(S3S0):
    double s3_s3s0 = simplex[3].dot(S0S3);
    //This dot product is equivalent to -S3.dot(S3S1):
    double s3_s3s1 = simplex[3].dot(S1S3);
    //This dot product is equivalent to -S3.dot(S3S2):
    double s3_s3s2 = simplex[3].dot(S2S3);

    //=============================
    //Is origin in Voroni region of veretex 3
    //hout << "Vertex 2: s3_s3s0=" << s3_s3s0 << " s3_s3s1=" << s3_s3s1 << " s3_s3s2=" << s3_s3s2 << endl;
    if (s3_s3s0 < Zero && s3_s3s1 < Zero && s3_s3s2 < Zero) {

        //Update simplex to cotain only vertex 3
        //hout << "Vertex 3 from tetrahedron" << endl;
        if (!Reduce_simplex_4_to_1(3, simplex)) {
            hout << "Error when reducing simplex from 4 to 1. Vertex 3." << endl;
            return -1;
        }

        //Terminate the function
        return 1;
    }

    //=============================
    //Edges
    //hout << "Edges + + + +" << endl;

    //Precompute some normal vectors to the faces of the tetrahedron
    //The last point in Normal_to_face_ouside_tetrahedron is not in
    //the face composed by the first three points
    Point_3D N013 = Normal_to_face_ouside_tetrahedron(simplex[0], simplex[1], simplex[3], simplex[2]);
    //hout << "N013=" << N013.str() << endl;
    Point_3D N012 = Normal_to_face_ouside_tetrahedron(simplex[0], simplex[1], simplex[2], simplex[3]);
    //hout << "N012=" << N012.str() << endl;

    //Vector normal to edge S0S1 on the plane of face S0S1S2
    Point_3D N01_012 = S0S1.cross(N012);
    //Make sure it goes in the direction opposite to vertex S2
    if (N01_012.dot(S0S2) > Zero) {
        //Normal vector goes towards S2, so reverse it
        N01_012 = N01_012 * (-1.0);
    }

    //Vector normal to edge S0S1 on the plane of face S0S1S3
    Point_3D N01_013 = S0S1.cross(N013);
    //Make sure it goes in the direction opposite to vertex S3
    if (N01_013.dot(S0S3) > Zero) {
        //Normal vector goes towards S3, so reverse it
        N01_013 = N01_013 * (-1.0);
    }

    //Precompute some dot products
    double s0_n01_012 = -simplex[0].dot(N01_012);
    double s0_n01_013 = -simplex[0].dot(N01_013);

    //=============================
    //Is origin in Voroni region of edge S0S1
    //hout << "Edge 01: s0_s0s1=" << s0_s0s1 << " s1_s1s0=" << s1_s1s0 << " s0_n01_012=" << s0_n01_012 << " s0_n01_013=" << s0_n01_013 << endl;
    if (s0_s0s1 > Zero && s1_s1s0 > Zero && s0_n01_012 > Zero && s0_n01_013 > Zero) {

        //Update simplex to contain only S0S1
        //hout << "Edge 01 from tetrahedron" << endl;
        if (!Reduce_simplex_4_to_2(0, 1, simplex)) {
            hout << "Error when reducing simplex from 4 to 2. Edge 01." << endl;
            return -1;
        }

        //Terminate the function
        return 1;
    }

    //Precompute some normal vectors to the faces of the tetrahedron
    //The last point in Normal_to_face_ouside_tetrahedron is not in
    //the face composed by the first three points
    Point_3D N023 = Normal_to_face_ouside_tetrahedron(simplex[0], simplex[2], simplex[3], simplex[1]);
    //hout << "N023=" << N023.str() << endl;

    //Vector normal to edge S0S2 on the plane of face S0S1S2
    Point_3D N02_012 = S0S2.cross(N012);
    //Make sure it goes in the direction opposite to vertex S1
    if (N02_012.dot(S0S1) > Zero) {
        //Normal vector goes towards S1, so reverse it
        N02_012 = N02_012 * (-1.0);
    }

    //Vector normal to edge S0S2 on the plane of face S0S2S3
    Point_3D N02_023 = S0S2.cross(N023);
    //Make sure it goes in the direction opposite to vertex S3
    if (N02_023.dot(S0S3) > Zero) {
        //Normal vector goes towards S3, so reverse it
        N02_023 = N02_023 * (-1.0);
    }

    //Precompute some dot products
    double s0_n02_012 = -simplex[0].dot(N02_012);
    double s0_n02_023 = -simplex[0].dot(N02_023);

    //=============================
    //Is origin in Voroni region of edge S0S2
    //hout << "Edge 02: s0_s0s2=" << s0_s0s2 << " s2_s2s0=" << s2_s2s0 << " s0_n02_012=" << s0_n02_012 << " s0_n02_023=" << s0_n02_023 << endl;
    if (s0_s0s2 > Zero && s2_s2s0 > Zero && s0_n02_012 > Zero && s0_n02_023 > Zero) {

        //Update simplex to contain only S0S2
        //hout << "Edge 02 from tetrahedron" << endl;
        if (!Reduce_simplex_4_to_2(0, 2, simplex)) {
            hout << "Error when reducing simplex from 4 to 2. Edge 02." << endl;
            return -1;
        }

        //Terminate the function
        return 1;
    }

    //Vector normal to edge S0S3 on the plane of face S0S1S3
    Point_3D N03_013 = S0S3.cross(N013);
    //Make sure it goes in the direction opposite to vertex S1
    if (N03_013.dot(S0S1) > Zero) {
        //Normal vector goes towards S1, so reverse it
        N03_013 = N03_013 * (-1.0);
    }

    //Vector normal to edge S0S3 on the plane of face S0S2S3
    Point_3D N03_023 = S0S3.cross(N023);
    //Make sure it goes in the direction opposite to vertex S2
    if (N03_023.dot(S0S2) > Zero) {
        //Normal vector goes towards S2, so reverse it
        N03_023 = N03_023 * (-1.0);
    }

    //Precompute some dot products
    double s0_n03_013 = -simplex[0].dot(N03_013);
    double s0_n03_023 = -simplex[0].dot(N03_023);

    //=============================
    //Is origin in Voroni region of edge S0S3
    //hout << "Edge 03: s0_s0s3=" << s0_s0s3 << " s3_s3s0=" << s3_s3s0 << " s0_n03_013=" << s0_n03_013 << " s0_n03_023=" << s0_n03_023 << endl;
    if (s0_s0s3 > Zero && s3_s3s0 > Zero && s0_n03_013 > Zero && s0_n03_023 > Zero) {

        //Update simplex to contain only S0S3
        //hout << "Edge 03 from tetrahedron" << endl;
        if (!Reduce_simplex_4_to_2(0, 3, simplex)) {
            hout << "Error when reducing simplex from 4 to 2. Edge 03." << endl;
            return -1;
        }

        //Terminate the function
        return 1;
    }

    //Precompute a normal vector to a face of the tetrahedron
    //The last point in Normal_to_face_ouside_tetrahedron is not in
    //the face composed by the first three points
    Point_3D N123 = Normal_to_face_ouside_tetrahedron(simplex[1], simplex[2], simplex[3], simplex[0]);
    //hout << "N123=" << N123.str() << endl;

    //Vector normal to edge S1S2 on the plane of face S0S1S2
    Point_3D N12_012 = S1S2.cross(N012);
    //Make sure it goes in the direction opposite to vertex S0
    if (N12_012.dot(S0S1) < Zero) { //This is equivalent to N12_012.dot(S1S0) > Zero
        //Normal vector goes towards S0, so reverse it
        N12_012 = N12_012 * (-1.0);
    }

    //Vector normal to edge S1S2 on the plane of face S1S2S3
    Point_3D N12_123 = S1S2.cross(N123);
    //Make sure it goes in the direction opposite to vertex S3
    if (N12_123.dot(S1S3) > Zero) {
        //Normal vector goes towards S3, so reverse it
        N12_123 = N12_123 * (-1.0);
    }

    //Precompute some dot products
    double s1_n12_012 = -simplex[1].dot(N12_012);
    double s1_n12_123 = -simplex[1].dot(N12_123);

    //=============================
    //Is origin in Voroni region of edge S1S2
    //hout << "Edge 12: s1_s1s2=" << s1_s1s2 << " s2_s2s1=" << s2_s2s1 << " s1_n12_012=" << s1_n12_012 << " s1_n12_123=" << s1_n12_123 << endl;
    if (s1_s1s2 > Zero && s2_s2s1 > Zero && s1_n12_012 > Zero && s1_n12_123 > Zero) {

        //Update simplex to contain only S1S2
        //hout << "Edge 12 from tetrahedron" << endl;
        if (!Reduce_simplex_4_to_2(1, 2, simplex)) {
            hout << "Error when reducing simplex from 4 to 2. Edge 12." << endl;
            return -1;
        }

        //Terminate the function
        return 1;
    }

    //Vector normal to edge S1S3 on the plane of face S0S1S3
    Point_3D N13_013 = S1S3.cross(N013);
    //Make sure it goes in the direction opposite to vertex S0
    if (N13_013.dot(S0S1) < Zero) { //This is equivalent to N13_013.dot(S1S0) > Zero
        //Normal vector goes towards S0, so reverse it
        N13_013 = N13_013 * (-1.0);
    }

    //Vector normal to edge S1S3 on the plane of face S1S2S3
    Point_3D N13_123 = S1S3.cross(N123);
    //Make sure it goes in the direction opposite to vertex S2
    if (N13_123.dot(S1S2) > Zero) {
        //Normal vector goes towards S2, so reverse it
        N13_123 = N13_123 * (-1.0);
    }

    //Precompute some dot products
    double s1_n13_013 = -simplex[1].dot(N13_013);
    double s1_n13_123 = -simplex[1].dot(N13_123);

    //=============================
    //Is origin in Voroni region of edge S1S3
    //hout << "Edge 13: s1_s1s3=" << s1_s1s3 << " s3_s3s1=" << s3_s3s1 << " s1_n13_013=" << s1_n13_013 << " s1_n13_123=" << s1_n13_123 << endl;
    if (s1_s1s3 > Zero && s3_s3s1 > Zero && s1_n13_013 > Zero && s1_n13_123 > Zero) {

        //Update simplex to contain only S1S3
        //hout << "Edge 13 from tetrahedron" << endl;
        if (!Reduce_simplex_4_to_2(1, 3, simplex)) {
            hout << "Error when reducing simplex from 4 to 2. Edge 13." << endl;
            return -1;
        }

        //Terminate the function
        return 1;
    }

    //Vector normal to edge S2S3 on the plane of face S0S2S3
    Point_3D N23_023 = S2S3.cross(N023);
    //Make sure it goes in the direction opposite to vertex S0
    if (N23_023.dot(S0S2) < Zero) { //This is equivalent to N23_023.dot(S2S0) > Zero
        //Normal vector goes towards S0, so reverse it
        N23_023 = N23_023 * (-1.0);
    }

    //Vector normal to edge S2S3 on the plane of face S1S2S3
    Point_3D N23_123 = S2S3.cross(N123);
    //Make sure it goes in the direction opposite to vertex S1
    if (N23_123.dot(S1S2) < Zero) { //This is equivalent to N23_123.dot(S2S1) > Zero
        //Normal vector goes towards S1, so reverse it
        N23_123 = N23_123 * (-1.0);
    }

    //Precompute some dot products
    double s2_n23_023 = -simplex[2].dot(N23_023);
    double s2_n23_123 = -simplex[2].dot(N23_123);

    //=============================
    //Is origin in Voroni region of edge S2S3
    //hout << "Edge 23: s2_s2s3=" << s2_s2s3 << " s3_s3s2=" << s3_s3s2 << " s2_n23_023=" << s2_n23_023 << " s2_n23_123=" << s2_n23_123 << endl;
    if (s2_s2s3 > Zero && s3_s3s2 > Zero && s2_n23_023 > Zero && s2_n23_123 > Zero) {

        //Update simplex to contain only S2S3
        //hout << "Edge 23 from tetrahedron" << endl;
        if (!Reduce_simplex_4_to_2(2, 3, simplex)) {
            hout << "Error when reducing simplex from 4 to 2. Edge 23." << endl;
            return -1;
        }

        //Terminate the function
        return 1;
    }

    //=============================
    //Faces
    //hout << "Faces + - - -" << endl;

    //=============================
    //Is origin in Voroni region of face S0S1S2
    //hout << "Face 012: -simplex[0].dot(N012)=" << -simplex[0].dot(N012);
    //hout << " s0_n01_012=" << s0_n01_012 << " s0_n02_012=" << s0_n02_012 << " s1_n12_012=" << s1_n12_012 << endl;
    if (-simplex[0].dot(N012) > Zero && s0_n01_012 < Zero && s0_n02_012 < Zero && s1_n12_012 < Zero) {

        //Update simplex to contain only S0S1S2
        //hout << "Face 012 from tetrahedron" << endl;
        simplex.pop_back();

        return 1;
    }

    //=============================
    //Is origin in Voroni region of face S0S2S3
    //hout << "Face 012: -simplex[0].dot(N023)=" << -simplex[0].dot(N023);
    //hout << " s0_n02_023=" << s0_n02_023 << " s0_n03_023=" << s0_n03_023 << " s2_n23_023=" << s2_n23_023 << endl;
    if (-simplex[0].dot(N023) > Zero && s0_n02_023 < Zero && s0_n03_023 < Zero && s2_n23_023 < Zero) {

        //Update simplex to contain only S0S2S3
        //hout << "Face 023 from tetrahedron" << endl;
        simplex[1] = simplex[3];
        simplex.pop_back();

        return 1;
    }

    //=============================
    //Is origin in Voroni region of face S0S1S3
    //hout << "Face 012: -simplex[0].dot(N013)=" << -simplex[0].dot(N013);
    //hout << " s0_n01_013=" << s0_n01_013 << " s1_n13_013=" << s1_n13_013 << " s0_n03_013=" << s0_n03_013 << endl;
    if (-simplex[0].dot(N013) > Zero && s0_n01_013 < Zero && s1_n13_013 < Zero && s0_n03_013 < Zero) {

        //Update simplex to contain only S0S1S3
        //hout << "Face 013 from tetrahedron" << endl;
        simplex[2] = simplex[3];
        simplex.pop_back();

        return 1;
    }

    //=============================
    //Is origin in Voroni region of face S1S2S3
    //hout << "Face 012: -simplex[1].dot(N123)=" << -simplex[1].dot(N123);
    //hout << " s1_n12_123=" << s1_n12_123 << " s1_n13_123=" << s1_n13_123 << " s2_n23_123=" << s2_n23_123 << endl;
    if (-simplex[1].dot(N123) > Zero && s1_n12_123 < Zero && s1_n13_123 < Zero && s2_n23_123 < Zero) {

        //Update simplex to contain only S1S2S3
        //hout << "Face 123 from tetrahedron" << endl;
        simplex[0] = simplex[3];
        simplex.pop_back();

        return 1;
    }

    //If this part of the code is reached, then the origin is inside the tetrahedron
    //This part of the code should not be reached, since this function should
    //only be used when it is known that there is no intersection of polyhedra
    //hout << "Origin is inside tetrahedron" << endl;

    return 0;
}
//This function removes all elements from a simplex except for one, which is indicated by an index
int Collision_detection::Reduce_simplex_4_to_1(const int& idx, vector<Point_3D>& simplex) 
{
    //Check for valid index
    if (idx < 0 || idx > 3) {
        hout << "Error in reduce_simplex_4_to_1. Invalid index. Index can be any integer in [0, 3]. Input was: " << idx << endl;
        return 0;
    }

    //Temporary point to store point in index idx
    Point_3D tmp = simplex[idx];

    //Empty simplex
    simplex.clear();

    //Add vertex stored in tmp
    simplex.push_back(tmp);

    return 1;
}
//This function removes all elements from a simplex except for two, which are indicated by indices
int Collision_detection::Reduce_simplex_4_to_2(const int& idx1, const int& idx2, vector<Point_3D>& simplex) 
{
    //Check for valid indices
    if (idx1 < 0 || idx1 > 3) {
        hout << "Error in reduce_simplex_4_to_2. Invalid index1. Index can be any integer in [0, 3]. Input was: " << idx1 << endl;
        return 0;
    }
    if (idx2 < 0 || idx2 > 3) {
        hout << "Error in reduce_simplex_4_to_2. Invalid index2. Index can be any integer in [0, 3]. Input was: " << idx2 << endl;
        return 0;
    }

    //Temporary points
    Point_3D tmp1 = simplex[idx1];
    Point_3D tmp2 = simplex[idx2];

    //Empty simplex
    simplex.clear();

    //Add stored vertices in temporary variables
    simplex.push_back(tmp1);
    simplex.push_back(tmp2);

    return 1;
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
    
    //hout<<"normal_and_distance_to_origin"<<endl;
    for (int i = 0; i < (int)Faces.size(); i++) {
        
        Normal_and_distance_to_origin(simplex, Faces[i], Normals[i], distances[i]);
        //hout<<"dist="<<distances[i]<<endl;
    }
    
    //Enter the loop to expand the polytope
    while (true) {
        
        //Find the face of the simplex closest to the origin
        //Actually I only need the normal and distance
        Point_3D N;
        double PD_;
        //hout<<"find_closest_face"<<endl;
        if (!Find_closest_face(Faces, Normals, distances, N, PD_))
        {
            hout << "Error in Extended Polytope Algorithm (EPA) when calling Find_closest_face." << endl;
            return 0;
        }
        
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
    //hout<<"r1="<<(simplex[f.v2] - simplex[f.v1]).str()<<" r2="<<(simplex[f.v3] - simplex[f.v1]).str()<<endl;
    //hout<<"N cross="<<N.str();
    //Make unit
    N.make_unit();
    //hout<<"N unit="<<N.str()<<endl;
    
    //Update the normal going from the origin towards the face
    //The vector v1 goes from the origin towards the face f,
    //then, the vector -v1 goes from the face towards the origin
    //If the dot product is positive, then both vectors "go towards the same direction"
    //So if the dot product of N and -v1 is positive, then N needs to be reversed
    //hout<<"dir check="<<dot(N, simplex[f.v1]*(-1.0))<<endl;
    //Also, save the value of the dot procduct as it is used to calculate the distance
    //from the origin to the face f
    double signed_dist = N.dot(simplex[f.v1]*(-1.0));
    normal = (signed_dist > Zero)? N*(-1.0) : N;
    
    //Calculate distance to origin
    //distance = abs(dot(normal, simplex[f.v1]));
    distance = abs(signed_dist);
    
    //hout<<"Face="<<f.str()<<" N="<<normal.str()<<" dist="<<distance<<endl;
}
//This function finds the face that is closest to the origin
//The output of this function is actually the normal vector and distance to origin of the closest face
int Collision_detection::Find_closest_face(const vector<TrFace> &Faces, const vector<Point_3D> &Normals, const vector<double> &distances, Point_3D &normal, double &PD) {

    //Check if there are any faces
    if (Faces.empty())
    {
        hout << "Error in Find_closest_face. No faces were found." << endl;
        hout << "Faces.size=" << Faces.size() << " Normals.size=" << Normals.size() << " distances.size=" << distances.size() << endl;
        return 0;
    }

    //Initialize the minimum distance (i.e., the penetration depth PD) with the first face
    //hout << "Faces.size=" << Faces.size() << " Normals.size=" << Normals.size() << " distances.size=" << distances.size() << endl;
    PD = distances[0];
    
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
    //hout<<"Closest Face="<<Faces[idx].str()<<endl;

    return 1;
}
//Reconstruct the polyhedron
void Collision_detection::Reconstruct_simplex(const vector<Point_3D> &simplex, const int &S, vector<TrFace> &Faces, vector<Point_3D> &Normals, vector<double> &distances) {
    
    //Vector to store the indices of the faces to remove
    vector<int> idxs;
    
    //Find the faces that are seen by the support point S, these face will be removed
    for (int i = 0; i < (int)Normals.size(); i++) {
        
        //The face is seen by the support point, S, if the vector that goes from the face to S
        //is in the same direction as the normal of the face
        //hout<<"Faces["<<i<<"]="<<Faces[i].str()<<" dot(S-v1,N)="<<dot(simplex[S]-simplex[Faces[i].v1], Normals[i])<<endl;
        if (Normals[i].dot(simplex[S]-simplex[Faces[i].v1]) > Zero ) {
            
            //Add the index of faces to remove
            idxs.push_back(i);
        }
    }
    //hout<<"idxs.size="<<idxs.size()<<endl;
    
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
        //hout<<"Deleted Face="<<Faces[idxs[i]].str()<<endl;
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
        //hout<<"New Face="<<Faces.back().str()<<" new_dist="<<distances.back()<<endl;
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
int Collision_detection::Distance_and_direction_from_simplex_to_origin(const vector<Point_3D> &simplex, Point_3D &N, double &dist)
{
    //Vector AO is used a few times
    Point_3D AO = simplex[0]*(-1);
    
    //Calculate the distance from simplex to origin depending on the type of simplex (vertex, edge, or face)
    if (simplex.size() == 3) {
        //Simplex is a face
        //hout << "simplex.size=3" << endl;
        
        //Get the normal vector
        N = (simplex[1] - simplex[0]).cross(simplex[2] - simplex[0]);
        //Make the vector unitary
        N.make_unit();
        
        //Make sure it goes towards the origin, e.g., in the same direction as AO
        double signed_dist = N.dot(simplex[0]);
        if (signed_dist < Zero) {
            N = N*(-1.0);
        }
        
        //Calculate distance to origin
        dist = abs(signed_dist);
    }
    else if (simplex.size() == 2) {
        //Simplex is an edge, where simplex = {A,B}
        //hout << "simplex.size=2" << endl;
        
        //AB = B - A
        Point_3D AB = simplex[1] - simplex[0];
        
        //dp^2 = (AO.dot(AB))^2/||AB||^2
        double dot_tmp = AO.dot(AB);
        //This is actually the projection squared, this is the quantity needed to obtain dist
        double dp_2 = dot_tmp*dot_tmp/AB.length2();
        
        //The distance of the edge to the origin is calculated using Pythagoras
        dist = sqrt(AO.length2() - dp_2);
        
        //The direction vector is the triple cross product
        N = (AB.cross(AO)).cross(AB);
        N.make_unit();
        
    }
    else if (simplex.size() == 1) {
        //Simplex is a vertex
        //hout << "simplex.size=1" << endl;
        
        //The direction vector goes from the vertex to origin, i.e., along AO
        N = AO.unit();
        
        //The distance is the length of the vector from origin to point
        dist = AO.length();
    }
    else {
        //This should not happen
        hout<<"Error in Distance_and_direction_from_simplex_to_origin. Simplex has an invalid number of elements: "<<simplex.size()<<" elements. Simplex can only have 1, 2, o 3 elements."<<endl;
        return 0;
    }
    return 1;
}
