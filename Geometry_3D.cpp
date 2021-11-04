//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Defeinition of geometry elements (point, line, plane) and operations with them
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Geometry_3D.h"

//---------------------------------------------------------------------------
//Constructor
Point_3D::Point_3D( double px, double py, double pz )
{
	x = px;
	y = py;
	z = pz;
}
Point_3D Point_3D::operator+( const Point_3D &pt )const
{
    Point_3D rp( x + pt.x, y + pt.y, z + pt.z );
    return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator-( const Point_3D &pt )const
{
    Point_3D rp( x - pt.x, y - pt.y, z - pt.z );
    return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator+( const double &d )const
{
    Point_3D rp( x + d, y + d, z + d );
    return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator-( const double &d )const
{
    Point_3D rp( x - d, y - d, z - d );
    return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator*( const double &d )const
{
    Point_3D rp( x*d, y*d, z*d );
    return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator/( const double &d )const
{
    Point_3D rp( x/d, y/d, z/d );
    return rp;
}
//---------------------------------------------------------------------------
bool Point_3D::operator==(const Point_3D &pt )const
{
	return (x==pt.x&&y==pt.y&&z==pt.z);
}
//---------------------------------------------------------------------------
bool Point_3D::operator!=( Point_3D &pt )const
{
	return (x!=pt.x||y!=pt.y||z!=pt.z);
}
//---------------------------------------------------------------------------
double Point_3D::distance_to(const Point_3D &pt )const
{
	double rv2 = (x-pt.x)*(x-pt.x)+(y-pt.y)*(y-pt.y)+(z-pt.z)*(z-pt.z);
	return sqrt(rv2);
}
//---------------------------------------------------------------------------
double Point_3D::distance_to(const double &px, const double &py, const double &pz)const
{
	double rv2 = (x-px)*(x-px)+(y-py)*(y-py)+(z-pz)*(z-pz);
	return sqrt(rv2);
}
//---------------------------------------------------------------------------
double Point_3D::squared_distance_to(const Point_3D &pt )const
{
    return (x-pt.x)*(x-pt.x)+(y-pt.y)*(y-pt.y)+(z-pt.z)*(z-pt.z);
}
//---------------------------------------------------------------------------
double Point_3D::squared_distance_to(const double &px, const double &py, const double &pz)const
{
    return (x-px)*(x-px)+(y-py)*(y-py)+(z-pz)*(z-pz);
}
//---------------------------------------------------------------------------
//Length of the vector from the origin to Point_3D
double Point_3D::length()const
{
    return sqrt(x*x + y*y + z*z);
}
//---------------------------------------------------------------------------
//Squared length of the vector from the origin to Point_3D
double Point_3D::length2()const
{
    return (x*x + y*y + z*z);
}
//---------------------------------------------------------------------------
Point_3D Point_3D::cross(const Point_3D &point)const
{
    double a = y*point.z - point.y*z;
    double b = -(x*point.z - point.x*z);
    double c = x*point.y - point.x*y;
    
    return Point_3D(a,b,c);
}
//---------------------------------------------------------------------------
double Point_3D::dot(const Point_3D &point)const
{
    return x*point.x + y*point.y + z*point.z;
}
//---------------------------------------------------------------------------
//Calculates the rotation of a point plus a displacement
Point_3D Point_3D::rotation(const MathMatrix &Matrix, const Point_3D &displacement)const
{
    //Rotated point
    //new_point = Matrix*point + displacement
    Point_3D new_point;
    
    //x-coordinate
    new_point.x = Matrix.element[0][0]*x + Matrix.element[0][1]*y + Matrix.element[0][2]*z + displacement.x;
    
    //y-coordinate
    new_point.y = Matrix.element[1][0]*x + Matrix.element[1][1]*y + Matrix.element[1][2]*z + displacement.y;
    
    //z-coordinate
    new_point.z = Matrix.element[2][0]*x + Matrix.element[2][1]*y + Matrix.element[2][2]*z + displacement.z;
    
    return new_point;

}
//---------------------------------------------------------------------------
//Return a unit vector in the direction of the origin to the point
Point_3D Point_3D::unit()
{
    double length = sqrt(x*x + y*y + z*z);
    //Create and return a unit vector
    return Point_3D(x/length, y/length, z/length);
}
//---------------------------------------------------------------------------
//Make the point a unit vector in the direction of the origin to the point
void Point_3D::make_unit()
{
    double length = sqrt(x*x + y*y + z*z);
    //Normalize components
    x = x/length;
    y = y/length;
    z = z/length;
}
//---------------------------------------------------------------------------
void Point_3D::set(const double &x_, const double &y_, const double &z_)
{
    //Set the components as given by the arguments of the function
    x = x_;
    y = y_;
    z = z_;
}
//---------------------------------------------------------------------------
void Point_3D::set(const Point_3D &P)
{
    //Set the components to be the same as those in point P
    x = P.x;
    y = P.y;
    z = P.z;
}
//---------------------------------------------------------------------------
//Make a string version of the point
string Point_3D::str()const
{
    string s = "(" + to_string(x) + ", " + to_string(y) + ", " + to_string(z) + ")";
    
    return s;
}
//---------------------------------------------------------------------------
//Make a string version of the point with the indicated precision
string Point_3D::str(const int &prec)const
{
    //Define a string stream
    stringstream ss;
    
    //Set the presicion of the string stream
    ss.precision(prec);
    
    //Send the components of the point to the sting stream
    //The std::fixed is needed so that all other double sent to the sting stream
    //also have the precision set above
    //Otherwise, without the std::fixed, all doubles that follow will have the default precision
    ss<<fixed<<x<<", "<<y<<", "<<z;
    
    return ss.str();
}
//---------------------------------------------------------------------------
//This function determines if a point is outside a cuboid
bool Point_3D::is_outside_cuboid(const struct cuboid &cub)const
{
    //Check if the point is outside the cuboid
    if ((x - cub.poi_min.x < Zero)||(x - cub.max_x > Zero)||(y - cub.poi_min.y < Zero)||(y - cub.max_y > Zero)||(z - cub.poi_min.z < Zero)||(z - cub.max_z > Zero))
        //The point is outside the cuboid, so return true
        return true;
    //If not outside the cuboid, then it is inside or at the boudary, so return false
    return false;
}
//---------------------------------------------------------------------------
//This function determines if a point is at a cuboid boundary
bool Point_3D::is_at_cuboid_boundary(const struct cuboid &cub)const
{
    //Check if the point is close enough to the cuboid boundary
    if ((abs(x - cub.poi_min.x) < Zero)||(abs(x - cub.max_x) < Zero)||(abs(y - cub.poi_min.y) < Zero)||(abs(y - cub.max_y) < Zero)||(abs(z - cub.poi_min.z) < Zero)||(abs(z - cub.max_z) < Zero))
        //The point is at the boundary
        return true;
    //If not at the cuboid boundary, then it is inside or outside the cuboid
    return false;
}
//---------------------------------------------------------------------------
//This function determines is a point is at a cuboid boundary
//===========================================================================
//Functions for the 3D Line class
//---------------------------------------------------------------------------
//Constructor
Line_3D::Line_3D(Point_3D p0, Point_3D p1)
{
	point[0] = p0;
	point[1] = p1;
	xm = p1.x - p0.x;
	yn = p1.y - p0.y;
	zl = p1.z - p0.z;
	len = length();
    if (len < Zero) virtual_line = false;
    else virtual_line = true;
}
//---------------------------------------------------------------------------   
//Get length of the line segment
double Line_3D::length() 
{
	double dx = point[1].x-point[0].x;
	double dy = point[1].y-point[0].y;
	double dz = point[1].z-point[0].z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}
//---------------------------------------------------------------------------
//Calculate the distance from the line segment to the specified point 
double Line_3D::distance_point_to_line(const Point_3D *point_temp)const
{
	double dis = 0;
    if (abs(xm) < Zero && abs(yn) < Zero && abs(zl) < Zero)
	{
		hout << "Warning: this segment is reduced to a point!" <<endl;
		double X = point_temp->x-point[0].x;
		double Y = point_temp->y-point[0].y;
		double Z = point_temp->z-point[0].z;
		dis = sqrt(X*X+Y*Y+Z*Z);
	}
	else
	{
		double X = yn*(point_temp->y-point[0].y)-zl*(point_temp->z-point[0].z);
		double Y = zl*(point_temp->z-point[0].z)-xm*(point_temp->x-point[0].x);
		double Z = xm*(point_temp->x-point[0].x)-yn*(point_temp->y-point[0].y);
		dis = sqrt(X*X+Y*Y+Z*Z)/sqrt(xm*xm+yn*yn+zl*zl);
	}
	return dis;
}
double Line_3D::distance_point_to_line(const Point_3D &point_temp)const
{
	double dis = 0;
    if (abs(xm) < Zero && abs(yn) < Zero && abs(zl) < Zero)
	{
		hout << "Warning: this segment is reduced to a point!" <<endl;
		double X = point_temp.x-point[0].x;
		double Y = point_temp.y-point[0].y;
		double Z = point_temp.z-point[0].z;
		dis = sqrt(X*X+Y*Y+Z*Z);
	}
	{
		double X = yn*(point_temp.y-point[0].y)-zl*(point_temp.z-point[0].z);
		double Y = zl*(point_temp.z-point[0].z)-xm*(point_temp.x-point[0].x);
		double Z = xm*(point_temp.x-point[0].x)-yn*(point_temp.y-point[0].y);
		dis = sqrt(X*X+Y*Y+Z*Z)/sqrt(xm*xm+yn*yn+zl*zl);
	}
	return dis;
}
double Line_3D::distance_point_to_line(const double dx, const double dy, const double dz)const
{
	double dis = 0;
    if (abs(xm) < Zero && abs(yn) < Zero && abs(zl) < Zero)
	{
		hout << "Warning: this segment is reduced to a point!" <<endl;
		double X = dx-point[0].x;
		double Y = dy-point[0].y;
		double Z = dz-point[0].z;
		dis = sqrt(X*X+Y*Y+Z*Z);
	}
	{
		double X = yn*(dy-point[0].y)-zl*(dz-point[0].z);
		double Y = zl*(dz-point[0].z)-xm*(dx-point[0].x);
		double Z = xm*(dx-point[0].x)-yn*(dy-point[0].y);
		dis = sqrt(X*X+Y*Y+Z*Z)/sqrt(xm*xm+yn*yn+zl*zl);
	}
	return dis;
}
//---------------------------------------------------------------------------
//Determine is a Point is on the line segment
int Line_3D::contain(const Point_3D &point_temp)const
{
	//to judge if the distance from a point to two endpoints is larger than the distance between endpoints
	if( fabs(point_temp.distance_to(point[0])+point_temp.distance_to(point[1])-point[0].distance_to(point[1]))>Zero ) return 0; 
	return 1;
}
//===========================================================================
//Functions for the 3D plane class
//---------------------------------------------------------------------------
//Constructor
Plane_3D::Plane_3D(double para[])
{
    //Copy the coefficients
	for(int i=0; i<4; i++)	coef[i] = para[i];
    
    //Calculate normal vector
    N.set(coef[0],coef[1],coef[2]);
    N.make_unit();
    
    if (abs(coef[0]) < Zero && abs(coef[1]) < Zero && abs(coef[2]) < Zero) virtual_plane = true;
	else virtual_plane = false;
}
//---------------------------------------------------------------------------
//Constructor
Plane_3D::Plane_3D(double a, double b, double c, double d)
{
    coef[0] = a;
    coef[1] = b;
    coef[2] = c;
    coef[3] = d;
    
    //Calculate normal vector
    N.set(coef[0],coef[1],coef[2]);
    N.make_unit();
    
    if (abs(coef[0]) < Zero && abs(coef[1]) < Zero && abs(coef[2]) < Zero) virtual_plane = true;
    else virtual_plane = false;
}
//---------------------------------------------------------------------------
//Constructor with three points
Plane_3D::Plane_3D(const Point_3D &P1, const Point_3D &P2, const Point_3D &P3)
{
    //Calculate normal vector and normalize it
    N = (P2 - P1).cross(P3 - P1);
    N.make_unit();
    
    //Top face
    coef[0] = N.x;
    coef[1] = N.y;
    coef[2] = N.z;
    coef[3] = N.dot(P1)*(-1);
    
    if (abs(coef[0]) < Zero && abs(coef[1]) < Zero && abs(coef[2]) < Zero) virtual_plane = true;
    else virtual_plane = false;
}
//---------------------------------------------------------------------------
//Determine if a point is on in this plane, where the point is given as a Point_3D object
int Plane_3D::contain(const Point_3D &point_temp)const
{
    if( abs(coef[0]*point_temp.x + coef[1]*point_temp.y + coef[2]*point_temp.z + coef[3]) < Zero) {
        //in the plane
        return 1;
    }
    //out of the plane
	return 0;
}
//---------------------------------------------------------------------------
//Determine if a point is on in this plane, where the point is given by its three components
int Plane_3D::contain(const double &dx, const double &dy, const double &dz)const
{
    if( abs(coef[0]*dx + coef[1]*dy + coef[2]*dz + coef[3]) < Zero ) {
        //in the plane
        return 1;
    }
    //out of the plane
	return 0;
}
//---------------------------------------------------------------------------
//Function to determine the distance from a plane to a point
double Plane_3D::distance_to(const Point_3D &P)const
{
    double num = coef[0]*P.x + coef[1]*P.y + coef[2]*P.z + coef[3];
    double den = sqrt(coef[0]*coef[0] + coef[1]*coef[1] + coef[2]*coef[2]);
    
    return abs(num/den);
}
//---------------------------------------------------------------------------
string Plane_3D::str() const
{
    string str = to_string(coef[0])+ "x ";
    if (coef[1] < Zero && abs(coef[1]) > Zero) {
        //Second coefficient is negative
        str += to_string(coef[1]);
    }
    else {
        //Second coefficient is positive
        str += "+ " + to_string(coef[1]);
    }
    str += "y ";
    if (coef[2] < Zero && abs(coef[2]) > Zero) {
        //Third coefficient is negative
        str += to_string(coef[2]);
    }
    else {
        //Third coefficient is positive
        str += "+ " + to_string(coef[2]);
    }
    str += "z ";
    if (coef[3] < Zero && abs(coef[3]) > Zero) {
        //Third coefficient is negative
        str += to_string(coef[3]);
    }
    else {
        //Third coefficient is positive
        str += "+ " + to_string(coef[3]);
    }
    str += " = 0";
    return str;
}
//===========================================================================
