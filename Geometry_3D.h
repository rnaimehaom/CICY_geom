//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Defeinition of geometry elements (point, line, plane) and operations with them
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef GEOMETRY_3D_H
#define GEOMETRY_3D_H

#include<cmath>
#include<stdlib.h>
#include<vector>
#include "Hns.h"
#include "MathMatrix.h"
using namespace hns;

const int MAX_INT = 65536; //2^16 for calculating a random number

//---------------------------------------------------------------------------
//Definition of 3D points
class Point_3D
{
public:
    //Point coordinates
    double x, y, z;
    //Flag, usually CNT or GNP number
    int flag;

    //Constructors
    Point_3D(){};
    Point_3D( double px, double py, double pz );

    //Operations with points
    Point_3D operator+( const Point_3D &pt )const;
    Point_3D operator-( const Point_3D &pt )const;
    Point_3D operator+( const double &d )const;
    Point_3D operator-( const double &d )const;
    Point_3D operator*( const double &d )const;
    Point_3D operator/( const double &d )const;
    bool operator==(const Point_3D &pt )const;
    bool operator!=( Point_3D &pt )const;
    double distance_to(const Point_3D &pt)const;
    double distance_to(const double &px, const double &py,  const double &pz)const;
    double squared_distance_to(const Point_3D &pt)const;
    double squared_distance_to(const double &px, const double &py,  const double &pz)const;
    double length()const;
    double length2()const;
    Point_3D cross(const Point_3D &point)const;
    double dot(const Point_3D &point)const;
    Point_3D rotation(const MathMatrix &Matrix, const Point_3D &displacement)const;
    Point_3D unit();
    void make_unit();
    void set(const double &x_, const double &y_, const double &z_);
    void set(const Point_3D &P);
    string str()const;
    string str(const int &prec)const;
    bool is_outside_cuboid(const struct cuboid &cub)const;
    bool is_at_cuboid_boundary(const struct cuboid &cub)const;
};
//---------------------------------------------------------------------------
//Definition for 3D line segements
class Line_3D
{
public:
    //
    //Coordinates of the segment's endpoints
    Point_3D point[2];
    //Segment length
    double len;
    //Flag to indicate if it is a virtual (false) segment
    //false: reduced to a point
    //true: a real segment
    bool virtual_line;

    //Constructor
    Line_3D(Point_3D p0, Point_3D p1);

    //Get segment length
    double length();
    //Distance from line segment to a Point (different functions depending how the point is defined)
    double distance_point_to_line(const Point_3D *point_temp)const;
    double distance_point_to_line(const Point_3D &point_temp)const;
    double distance_point_to_line(const double dx, const double dy, const double dz)const;
    //Determine is a Point is on the line segment
    int contain(const Point_3D &point_temp)const;

private:
    //Coefficients from the equation of a line such that: (x-x0)/xm=(y-y0)/yn=(z-z0)/zl
    double xm, yn, zl ;
};
//---------------------------------------------------------------------------
//Definition for a plane in 3D
class Plane_3D
{
public:

    //Four coefficients that define the equation of plane: ax+by+cz+d=0
    //coef[0] = a, coef[1] = b, etc.
    double coef[4];
    //Normal unit vector to the plane
    Point_3D N;
    
    //flag to define if it is a virtual (false) plane
    //false: for virtual plane, its normal vector is (0,0,0)
    //true: for real plane
    bool virtual_plane;

    //Constructor
    Plane_3D(){};
    Plane_3D(double para[]);
    Plane_3D(double a, double b, double c, double d);
    Plane_3D(const Point_3D &P1, const Point_3D &P2, const Point_3D &P3);
    
    //Determine if a point is on this plane, where the point is given as a Point_3D object
    int contain(const Point_3D &point_temp)const;
    //Determine if a point is on this plane, where the point is given by its three components
    int contain(const double &dx, const double &dy, const double &dz)const;
    double distance_to(const Point_3D &P)const;
    string str() const;
};
//---------------------------------------------------------------------------
//Data structure for a cuboid
struct cuboid
{
    //"Origin" point
    Point_3D poi_min;
    //Cuboid's length, width and height
    double len_x, wid_y, hei_z;
    //Cuboid's maximum corrdinates
    double max_x, max_y, max_z;
    
    string str() const
    {
        //Define a string stream
        stringstream ss;
        
        //Add elements to string stream
        ss<<"P="<<poi_min.str()<<" max_x="<<max_x<<" max_y="<<max_y<<" max_z="<<max_z<<endl;
        
        return ss.str();
    }
};
//---------------------------------------------------------------------------
//Data structure for a junction
struct Junction
{
    //Point number on particle 1
    long int P1;
    //Type of particle 1, "CNT" or "GNP"
    string type1;
    //Particle number
    int N1;
    //Point number on particle 2
    long int P2;
    //Type of particle 1, "CNT" or "GNP"
    string type2;
    //Particle number
    int N2;
    //Junction distance
    double junction_dist;
    
    Junction(){}
    Junction(const long int &p1, const int &n1, const string &str1, const long int &p2, const int &n2, const string &str2, const double &jd){
        P1 = p1;
        P2 = p2;
        N1 = n1;
        N2 = n2;
        type1 = str1;
        type2 = str2;
        junction_dist = jd;
    }
};
//---------------------------------------------------------------------------
//Data structure for a shel
struct Shell {
    int shell_min;
    int shell_max;
};
//---------------------------------------------------------------------------
//Data structure for a triangular face based on the indices of a vector of Point_3D objects
struct TrFace {
    
    //Indices of the three vertices of the triangular face
    int v1, v2, v3;
    
    //Default contructor
    TrFace () {
        v1 = 0;
        v2 = 0;
        v3 = 0;
    }
    TrFace (const int &v1_, const int &v2_, const int &v3_) {
        v1=v1_;v2=v2_;v3=v3_;
    }
    
    void set(const int &v1_, const int &v2_, const int &v3_) {
        v1=v1_;v2=v2_;v3=v3_;
    }
    void set(const TrFace &f) {
        v1=f.v1;v2=f.v2;v3=f.v3;
    }
    string str(){
        return ("("+to_string(v1)+", "+to_string(v2)+", "+to_string(v3)+")");
    }
    string str()const{
        return ("("+to_string(v1)+", "+to_string(v2)+", "+to_string(v3)+")");
    }
};
//---------------------------------------------------------------------------
//Data structure for a triangular face based on the indices of a vector of Point_3D objects
//This TrFace data structure uses long ints to store the vertex numbers
struct TrFaceL {
    
    //Indices of the three vertices of the triangular face
    long int v1, v2, v3;
    
    //Default contructor
    TrFaceL () {
        v1 = 0;
        v2 = 0;
        v3 = 0;
    }
    TrFaceL (const long int &v1_, const long int &v2_, const long int &v3_) {
        v1=v1_;v2=v2_;v3=v3_;
    }
    
    void set(const long int &v1_, const long int &v2_, const long int &v3_) {
        v1=v1_;v2=v2_;v3=v3_;
    }
    void set(const TrFaceL &f) {
        v1=f.v1;v2=f.v2;v3=f.v3;
    }
    string str(){
        return ("("+to_string(v1)+", "+to_string(v2)+", "+to_string(v3)+")");
    }
    string str()const{
        return ("("+to_string(v1)+", "+to_string(v2)+", "+to_string(v3)+")");
    }
};
//---------------------------------------------------------------------------
//Data structure for an edge based on the indices of a vector of Point_3D objects
struct Edge {
    
    //Vertices of the edge
    int v1, v2;
    
    //Default contructor
    Edge () {
        v1 = 0;
        v2 = 0;
    }
    Edge (const int &v1_, const int &v2_) {
        v1=v1_;v2=v2_;
    }
    
    void set(const int &v1_, const int &v2_) {
        v1=v1_;v2=v2_;
    }
    void set(const Edge &e) {
        v1=e.v1;v2=e.v2;
    }
    string str(){
        return ("("+to_string(v1)+", "+to_string(v2)+")");
    }
    string str()const{
        return ("("+to_string(v1)+", "+to_string(v2)+")");
    }
    //Comparing two edges, here edges are not directed so AB and BA are the same edge
    bool operator==(const Edge &e) const {
        return ( (e.v1 == v1 && e.v2 == v2) || (e.v2 == v1 && e.v1 == v2) );
    }
};
//---------------------------------------------------------------------------
//Data structure for an edge based on the indices of a vector of Point_3D objects
//This Edge data structure uses long ints to store the vertex numbers
struct EdgeL {
    
    //Vertices of the edge
    long int v1, v2;
    
    //Default contructor
    EdgeL () {
        v1 = 0;
        v2 = 0;
    }
    EdgeL (const long int &v1_, const long int &v2_) {
        v1=v1_;v2=v2_;
    }
    
    void set(const long int &v1_, const long int &v2_) {
        v1=v1_;v2=v2_;
    }
    void set(const EdgeL &e) {
        v1=e.v1;v2=e.v2;
    }
    string str(){
        return (to_string(v1)+"-"+to_string(v2));
    }
    string str()const{
        return (to_string(v1) + "-" + to_string(v2));
    }
    //Comparing two edges, here edges are not directed so AB and BA are the same edge
    bool operator==(const EdgeL &e) const {
        return ( (e.v1 == v1 && e.v2 == v2) || (e.v2 == v1 && e.v1 == v2) );
    }
};
//---------------------------------------------------------------------------
//Data structure for GNPs
struct GNP {
    
    //Vertices of the GNP
    Point_3D vertices[8];
    //Plane equations for the GNP faces
    Plane_3D faces[6];
    //GNP center
    Point_3D center;
    //Rotation matrix
    MathMatrix rotation;
    //Side lengths
    double l;
    //Thickness
    double t;
    //Volume
    double volume;
    //Tringulation edges
    vector<EdgeL> triangulation;
    //Flag for GNP number
    int flag;
    
    //Default constructor
    GNP() {
        
        //Initialize rotation matrix
        rotation = MathMatrix(3,3);
    }
    
    //Function to determine if a point is inside the GNP
    bool Is_point_inside_gnp(const Point_3D &P)const
    {
        //This vector is used multiple times
        Point_3D V4P = P - vertices[4];
        
        //Check if P is between faces 0 and 1
        if (faces[0].N.dot(P - vertices[0]) < Zero && faces[1].N.dot(V4P) < Zero) {
            
            //Check if P is between faces 2 and 4
            if (faces[2].N.dot(V4P) < Zero && faces[4].N.dot(P - vertices[1]) < Zero) {
                
                //Check if P is between faces 3 and 5
                if (faces[3].N.dot(V4P) < Zero && faces[5].N.dot(P - vertices[2]) < Zero) {
                    
                    //Point P is inside the GNP so return true
                    return true;
                }
            }
        }
        
        //If this part of the code is reached, then the point P is outside the GNP
        return false;
    }
    
};

//---------------------------------------------------------------------------
//Data structure for Hybrid particles
struct Hybrid {
    
    //GNP number
    int GNP;
    //CNTs attached to the GNP
    vector<int> cnts_top, cnts_bottom;
    
    //Default constructor
    Hybrid() {}
};
#endif
//===========================================================================
