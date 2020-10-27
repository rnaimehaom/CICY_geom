//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Export visualization files in ASCII VTK legacy file format
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "VTK_Export.h"

//This function exports a VTK visualization file for CNTs when the points are in a 1D vector
int VTK_Export::Export_cnts_1D_vector(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const string &filename)const
{
    //Check that the points vector has points
    if (points.empty()) {
        hout<<"Vector of points is empty. NO visualization file was exported."<<endl;
        return 1;
    }
    
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_cnts_1D_vector when calling Add_header"<<endl;
        return 0;
    }
    
    //Add the line with the number of points
    otec<<"POINTS "<<points.size()<<" float" <<endl;
    
    //Add all the points
    if (!Add_points_from_vector(points, otec)) {
        hout<<"Error in Export_cnts_1D_vector when calling Add_points_from_vector"<<endl;
        return 0;
    }
    
    //Add the line indicating the number of lines+1, and the number of points in those lines
    otec<<"LINES "<<structure.size()+1<<' '<<points.size()<<endl;
    
    //Add the offsets:
    //The number of points used after adding each line, starting with a zero
    if (!Add_offsets_for_structure(structure, otec)) {
        hout<<"Error in Export_cnts_1D_vector when calling Add_offsets_for_structure"<<endl;
        return 0;
    }
    
    //Add connectivity
    if (!Add_connectivity_for_structure(structure, otec)) {
        hout<<"Error in Export_cnts_1D_vector when calling Add_connectivity_for_structure"<<endl;
        return 0;
    }
    
    //Close the file
    otec.close();
    
    return 1;
}
//This function adds the first four lines of the VTK file, which are always the same
//independently of the type of nanoparticle that needs to be exported
int VTK_Export::Add_header(ofstream &otec)const
{
    otec << "# vtk DataFile Version 5.1"<<endl;
    otec << "vtk output"<<endl;
    otec << "ASCII"<<endl;
    otec << "DATASET POLYDATA"<<endl;
    
    return 1;
}
//This function adds the coordinates of a vector of points into a file
//Four points are printed per line
int VTK_Export::Add_points_from_vector(const vector<Point_3D> &points, ofstream &otec)const
{
    
    //Add the point coordinates, separated by spaces and in groups of four points
    //(12 coordinates) per line
    
    //Add the first point
    otec<<points[0].x<<' '<<points[0].y<<' '<<points[0].z<<' ';
    
    //Add the remaninig points
    for (size_t i = 1; i < points.size(); i++) {
        
        //Check if a new line needs to be started
        if (!(i%4)) {
            otec<<endl;
        }
        
        //Add point i
        otec<<points[i].x<<' '<<points[i].y<<' '<<points[i].z<<' ';
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
//This function adds the offsets when the structure vector is available
int VTK_Export::Add_offsets_for_structure(const vector<vector<long int> > &structure, ofstream &otec)const
{
    //Add the line with the OFFSETS command
    otec<<"OFFSETS vtktypeint64"<<endl;
    
    //Output a zero, which is required
    otec<<"0 ";
    
    //Variable to store the accumulated number of points, initialize with the first line
    long int acc_points = (long int)structure[0].size();
    
    //Output the accumulated number of points
    otec<<acc_points<<' ';
    
    //Output the accumulated number of points for the remaining lines
    //Print 20 per line
    for (size_t i = 1; i < structure.size(); i++) {
        
        //Check if a new line needs to be started
        if (!(i%20)) {
            otec<<endl;
        }
        
        //Add the number of points from line i to the accumulated number of points
        acc_points += (long int)structure[i].size();
        
        //Output the accumulated number of points
        otec<<acc_points<<' ';
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
//This function adds the connectivity when the structure vector is available
int VTK_Export::Add_connectivity_for_structure(const vector<vector<long int> > &structure, ofstream &otec)const
{
    //Add the line with the CONNECTIVITY command
    otec<<"CONNECTIVITY vtktypeint64"<<endl;
    
    //Variable to count the number of points
    long int n_points = 0;
    
    //Add the connectivity, one connectivity per line
    for (size_t i = 0; i < structure.size(); i++) {
        
        //Add the connectivity of line i
        for (size_t j = 0; j < structure[i].size(); j++) {
            
            //Add the consecutive number of point j in line i
            otec<<n_points<<' ';
            
            //Increase the count of points
            n_points++;
        }
        //Start a new line
        otec<<endl;
    }
    
    return 1;
}
//This function exports a VTK visualization file for CNTs when the points are in a 2D vector
int VTK_Export::Export_cnts_2D_vector(const vector<vector<Point_3D> > &points, const string &filename)const
{
    //Check that the points vector has points
    if (points.empty()) {
        hout<<"Vector of points is empty. NO visualization file was exported."<<endl;
        return 1;
    }
    
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_cnts_2D_vector when calling Add_header"<<endl;
        return 0;
    }
    
    //Count the number of points
    long int n_points = Count_number_of_points(points);
    
    //Add the line with the number of points
    otec<<"POINTS "<<n_points<<" float" <<endl;
    
    //Add all the points
    if (!Add_points_from_2D_vector(points, otec)) {
        hout<<"Error in Export_cnts_2D_vector when calling Add_points_from_2D_vector"<<endl;
        return 0;
    }
    
    //Add the line indicating the number of lines+1, and the number of points in those lines
    otec<<"LINES "<<points.size()+1<<' '<<n_points<<endl;
    
    //Add the offsets:
    //The number of points used after adding each line, starting with a zero
    if (!Add_offsets_for_2D_vector(points, otec)) {
        hout<<"Error in Export_cnts_2D_vector when calling Add_offsets_for_2D_vector"<<endl;
        return 0;
    }
    
    //Add connectivity
    if (!Add_connectivity_for_2D_vector(points, otec)) {
        hout<<"Error in Export_cnts_2D_vector when calling Add_connectivity_for_2D_vector"<<endl;
        return 0;
    }
    
    //Close the file
    otec.close();
    
    return 1;
}
//This function counts the number of points in a 2D vector
long int VTK_Export::Count_number_of_points(const vector<vector<Point_3D> > &points)const
{
    //Variable to store the number of points
    long int n_points = 0;
    
    //Count the number of points
    for (size_t i = 0; i < points.size(); i++) {
        
        //Add the number of points of CNT i to the total
        n_points += (long int)points[i].size();
    }
    
    return n_points;
}
//This function adds the coordinates of a 2D vector of points into a file
//For each vector of points, four points are printed per line
int VTK_Export::Add_points_from_2D_vector(const vector<vector<Point_3D> > &points, ofstream &otec)const
{
    
    //Add points, CNT by CNT
    for (size_t i = 0; i < points.size(); i++) {
        
        //Add points using the function for 1D vector
        if (!Add_points_from_vector(points[i], otec)) {
            hout<<"Error in Add_points_from_2D_vector when calling Add_points_from_vector"<<endl;
            return 0;;
        }
    }
    
    return 1;
}
//This function adds the offsets when the 2D vector of points is available
int VTK_Export::Add_offsets_for_2D_vector(const vector<vector<Point_3D> > &points, ofstream &otec)const
{
    //Add the line with the OFFSETS command
    otec<<"OFFSETS vtktypeint64"<<endl;
    
    //Output a zero, which is required
    otec<<"0 ";
    
    //Variable to store the accumulated number of points, initialize with the first line
    long int acc_points = (long int)points[0].size();
    
    //Output the accumulated number of points
    otec<<acc_points<<' ';
    
    //Output the accumulated number of points for the remaining lines
    //Print 20 per line
    for (size_t i = 1; i < points.size(); i++) {
        
        //Check if a new line needs to be started
        if (!(i%20)) {
            otec<<endl;
        }
        
        //Add the number of points from line i to the accumulated number of points
        acc_points += (long int)points[i].size();
        
        //Output the accumulated number of points
        otec<<acc_points<<' ';
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
//This function adds the connectivity when the 2D vector of points is available
int VTK_Export::Add_connectivity_for_2D_vector(const vector<vector<Point_3D> > &points, ofstream &otec)const
{
    //Add the line with the CONNECTIVITY command
    otec<<"CONNECTIVITY vtktypeint64"<<endl;
    
    //Variable to count the number of points
    long int n_points = 0;
    
    //Add the connectivity, one connectivity per line
    for (size_t i = 0; i < points.size(); i++) {
        
        //Add the connectivity of line i
        for (size_t j = 0; j < points[i].size(); j++) {
            
            //Add the consecutive number of point j in line i
            otec<<n_points<<' ';
            
            //Increase the count of points
            n_points++;
        }
        //Start a new line
        otec<<endl;
    }
    
    return 1;
}
//This function exports a vector of GNPs
int VTK_Export::Export_gnps(const vector<GNP> &gnps, const string &filename)const
{
    //Check that the points vector has points
    if (gnps.empty()) {
        hout<<"Vector of GNPs is empty. NO visualization file was exported."<<endl;
        return 1;
    }
    
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_gnps when calling Add_header"<<endl;
        return 0;
    }
    
    //Add the line with the number of points, which is the number of GNPs multiplied by 8
    otec<<"POINTS "<<gnps.size()*8<<" float" <<endl;
    
    //Add all the points, i.e., all the vertices of the GNPs
    if (!Add_all_gnp_vertices(gnps, otec)) {
        hout<<"Error in Export_gnps when calling Add_all_gnp_vertices"<<endl;
        return 0;
    }
    
    //Add the line with the polygons command, with the number of faces+1 and
    //the number of points in those faces
    otec<<"POLYGONS "<<gnps.size()*6 + 1<<' '<<gnps.size()*24<<endl;
    
    //Add the offsets
    if (!Add_ofsets_for_gnps(gnps, otec)) {
        hout<<"Error in Export_gnps when calling Add_ofsets_for_gnps"<<endl;
        return 0;
    }
    
    //Add the connectivity
    if (!Add_connectivity_for_gnps(gnps, otec)) {
        hout<<"Error in Export_gnps when calling Add_connectivity_for_gnps"<<endl;
        return 0;
    }
    
    //Close the file
    otec.close();
    
    return 1;
}
//This function adds all GNP vertices from all GNPs in a GNP vector
int VTK_Export::Add_all_gnp_vertices(const vector<GNP> &gnps, ofstream &otec)const
{
    
    //Add the points in the array of vertices for each GNP
    for (size_t i = 0; i < gnps.size(); i++) {
        
        //Add the vertices
        if (!Add_points_from_array(gnps[i].vertices, 8, otec)) {
            hout<<"Error in Add_all_gnp_vertices when calling Add_points_from_array"<<endl;
            return 0;
        }
        
        //Add a new line
        otec<<endl;
    }
    
    return 1;
}
//This function adds the coordinates of an array of points into a file
//Four points are printed per line
int VTK_Export::Add_points_from_array(const Point_3D points[], const int &arr_size, ofstream &otec)const
{
    //Add the point coordinates, separated by spaces and in groups of four points
    //(12 coordinates) per line
    
    //Add the first point
    otec<<points[0].x<<' '<<points[0].y<<' '<<points[0].z<<' ';
    
    //Add the remaninig points
    for (int i = 1; i < arr_size; i++) {
        
        //Check if a new line needs to be started
        if (!(i%4)) {
            otec<<endl;
        }
        
        //Add point i
        otec<<points[i].x<<' '<<points[i].y<<' '<<points[i].z<<' ';
    }
    
    return 1;
}
int VTK_Export::Add_ofsets_for_gnps(const vector<GNP> &gnps, ofstream &otec)const
{
    //Add the line with the OFFSETS command
    otec<<"OFFSETS vtktypeint64"<<endl;
    
    //Output a zero, which is required
    otec<<"0 ";
    
    //Variable to count the number of vertices
    int n_vertices = 0;
    
    //Write the offsets, one line per GNP
    for (size_t i = 0; i < gnps.size(); i++) {
        
        //Add the number of vertices per face in GNP i
        for (int j = 0; j < 6; j++) {
            
            //Add the four vertices of face j
            n_vertices += 4;
            otec<<n_vertices<<' ';
        }
        
        //Start a new line
        otec<<endl;
    }
    
    return 1;
}
int VTK_Export::Add_connectivity_for_gnps(const vector<GNP> &gnps, ofstream &otec, const long int &cnt_points_offset)const
{
    //Add the line with the CONNECTIVITY command
    otec<<"CONNECTIVITY vtktypeint64"<<endl;
    
    //Array with the indices of the vertices for each of the four faces
    long int faces[6][4] = {
        {0,1,2,3},
        {4,5,6,7},
        {3,0,4,7},
        {0,1,5,4},
        {1,2,6,5},
        {2,3,7,6},
    };
    
    for (long int i = 0; i < (int)gnps.size(); i++) {
        
        //Add each face j of GNP i
        for (int j = 0; j < 6; j++) {
            
            //Add vertex k of face j of GNP i
            for (int k = 0; k < 4; k++) {
                otec<<faces[j][k]+8*i+cnt_points_offset<<' ';
            }
        }
        
        //Start a new line
        otec<<endl;
    }
    
    return 1;
}
int VTK_Export::Export_hybrid_material(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const vector<GNP> &gnps, const string &filename)const
{
    //Check that at least one vector has elements
    if (points.empty() || gnps.empty()) {
        hout<<"Either the points vector (size = "<<points.size()<<") or the GNPs vector (size = "<<gnps.size()<<") is empty. NO visualization file was exported."<<endl;
        return 1;
    }
    
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_hybrid_material when calling Add_header"<<endl;
        return 0;
    }
    
    //Add the line with the number of points
    otec<<"POINTS "<<points.size() + gnps.size()*8<<" float" <<endl;
    
    //Add all CNT points first
    if (!Add_points_from_vector(points, otec)) {
        hout<<"Error in Export_hybrid_material when calling Add_points_from_vector"<<endl;
        return 0;
    }
    
    //Add all GNP points
    if (!Add_all_gnp_vertices(gnps, otec)) {
        hout<<"Error in Export_hybrid_material when calling Add_all_gnp_vertices"<<endl;
        return 0;
    }
    
    //---------------------------------------------------------------------------------------
    //Add the lines for CNTs
    
    //Add the line indicating the number of lines+1, and the number of points in those lines
    otec<<"LINES "<<structure.size()+1<<' '<<points.size()<<endl;
    
    //Add the offsets:
    //The number of points used after adding each line, starting with a zero
    if (!Add_offsets_for_structure(structure, otec)) {
        hout<<"Error in Export_hybrid_material when calling Add_offsets_for_structure"<<endl;
        return 0;
    }
    
    //Add connectivity
    if (!Add_connectivity_for_structure(structure, otec)) {
        hout<<"Error in Export_hybrid_material when calling Add_connectivity_for_structure"<<endl;
        return 0;
    }
    
    //---------------------------------------------------------------------------------------
    //Add the polygons for GNPs
    
    //Add the line with the polygons command, with the number of faces+1 and
    //the number of points in those faces
    otec<<"POLYGONS "<<gnps.size()*6 + 1<<' '<<gnps.size()*24<<endl;
    
    //Add the offsets
    if (!Add_ofsets_for_gnps(gnps, otec)) {
        hout<<"Error in Export_hybrid_material when calling Add_ofsets_for_gnps"<<endl;
        return 0;
    }
    
    //Add the connectivity
    //the first GNP point number is points.size(), i.e., the offset is points.size()
    if (!Add_connectivity_for_gnps(gnps, otec, (long int)points.size())) {
        hout<<"Error in Export_hybrid_material when calling Add_connectivity_for_gnps"<<endl;
        return 0;
    }
    
    return 1;
}
