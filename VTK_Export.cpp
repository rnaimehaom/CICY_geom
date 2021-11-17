//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Export visualization files in ASCII VTK legacy file format
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "VTK_Export.h"

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//CNTs
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
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
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int VTK_Export::Export_from_cnt_indices(const vector<Point_3D> &points, const vector<vector<long int> > &indices, const string &filename)const
{
    //Variable to store the number of points
    long int n_points = 0;
    //Variable to store the number of lines (CNTs)
    int n_cnts = 0;
    
    //Count the number of points and lines (CNTs)
    if (!Count_points_and_lines_from_indices(indices, n_points, n_cnts)) {
        hout<<"Error in Export_from_cnt_indices when calling Count_points_and_lines_from_indices"<<endl;
        return 0;
    }
    
    //Check that the points vector has points
    if (!n_points) {
        hout<<"No points to export. NO visualization file was exported for CNTs."<<endl;
        return 1;
    }
    
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_from_cnt_indices when calling Add_header"<<endl;
        return 0;
    }
    
    //Add the line with the number of points
    otec<<"POINTS "<<n_points<<" float" <<endl;
    
    //Add all the points
    if (!Add_points_from_indices(points, indices, otec)) {
        hout<<"Error in Export_from_cnt_indices when calling Add_points_from_indices"<<endl;
        return 0;
    }
    
    //Add the line indicating the number of lines+1, and the number of points in those lines
    otec<<"LINES "<<n_cnts+1<<' '<<n_points<<endl;
    
    //Add the offsets:
    //The number of points used after adding each line, starting with a zero
    if (!Add_offsets_from_indices(indices, otec)) {
        hout<<"Error in Export_from_cnt_indices when calling Add_offsets_from_indices"<<endl;
        return 0;
    }
    
    //Add connectivity
    if (!Add_connectivity_from_indices(indices, otec)) {
        hout<<"Error in Export_from_cnt_indices when calling Add_connectivity_from_indices"<<endl;
        return 0;
    }
    
    //Close the file
    otec.close();
    
    return 1;
}
int VTK_Export::Count_points_and_lines_from_indices(const vector<vector<long int> > &indices, long int &n_points, int &n_lines)const
{
    //Iterate over all the pairs of indices
    for (int i = 0; i < (int)indices.size(); i++) {
        for (int j = 0; j < (int)indices[i].size(); j=j+2) {
            
            //Add the difference of indices plus one
            n_points = n_points + indices[i][j+1] - indices[i][j] + 1;
            
            //Each pair of indices is a line (CNT), so increase the number of lines
            //each time j is increased by 2
            n_lines++;
        }
    }
    
    return 1;
}
//This function adds the coordinates of a vector of points into a file as indicated
//by a vector of indices
//Four points are printed per line
int VTK_Export::Add_points_from_indices(const vector<Point_3D> &points, const vector<vector<long int> > &indices, ofstream &otec)const
{
    //Add the point coordinates, separated by spaces and in groups of four points
    //(12 coordinates) per line
    
    //Variable to count the points
    long int n_points = 1;
    
    //Iterate over all the pairs of indices
    for (int i = 0; i < (int)indices.size(); i++) {
        for (int j = 0; j < (int)indices[i].size(); j= j+2) {
            
            //Check if a new line needs to be started
            if (!(n_points%4)) {
                otec<<endl;
            }
            
            //Get initial and final points of the line (CNT)
            long int P1 = indices[i][j];
            long int P2 = indices[i][j+1];
            
            //Iterate over all the points of the line (CNT segment) indicated
            //by the pair of indices j,j+1
            for (long int k = P1; k <= P2; k++) {
                
                //Append point i to file
                otec<<points[k].x<<' '<<points[k].y<<' '<<points[k].z<<' ';
                
                //Increase the count of points (for adding a new line only)
                n_points++;
            }
        }
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
int VTK_Export::Add_offsets_from_indices(const vector<vector<long int> > &indices, ofstream &otec)const
{
    //Add the line with the OFFSETS command
    otec<<"OFFSETS vtktypeint64"<<endl;
    
    //Output a zero, which is required
    otec<<"0 ";
    
    //Accumulate the number of points
    long int acc_points = 0;
    
    //Count the number of offsets
    int n_offsets = 1;
    
    //Iterate over all the pairs of indices
    for (int i = 0; i < (int)indices.size(); i++) {
        for (int j = 0; j < (int)indices[i].size(); j=j+2) {
            
            //Check if a new line needs to be started
            if (!(n_offsets%20)) {
                otec<<endl;
            }
            
            //Add the difference of indices plus one
            acc_points = acc_points + indices[i][j+1] - indices[i][j] + 1;
            
            //Output the accumulated number of points
            otec<<acc_points<<' ';
            
            //Each pair of indices is a line (CNT), so increase the number of lines
            //each time j is increased by 2
            n_offsets++;
        }
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
int VTK_Export::Add_connectivity_from_indices(const vector<vector<long int> > &indices, ofstream &otec)const
{
    //Add the line with the CONNECTIVITY command
    otec<<"CONNECTIVITY vtktypeint64"<<endl;
    
    //Variable to count the number of points
    long int n_points = 0;
    
    //Add the connectivity, one connectivity per line
    for (int i = 0; i < (int)indices.size(); i++) {
        for (int j = 0; j < (int)indices[i].size(); j= j+2) {
            
            //Get initial and final points of the line (CNT)
            long int P1 = indices[i][j];
            long int P2 = indices[i][j+1];
            
            //Iterate over all the points of the line (CNT segment) indicated
            //by the pair of indices j,j+1
            for (long int k = P1; k <= P2; k++) {
                
                //Add the consecutive number of point k in a line
                otec<<n_points<<' ';
                
                //Increase the number of points
                n_points++;
            }
            
            //For every pair of indices (i.e., for every iteration over j) start a new line
            otec<<endl;
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int VTK_Export::Export_cnts_in_cluster(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const vector<int> &cluster, const string &filename)const
{
    //Check that there are CNTs to export
    if (cluster.empty()) {
        hout<<"No CNTs to export. NO visualization file was exported for CNTs."<<endl;
        return 1;
    }
    
    //Variable to store the number of points
    long int n_points = 0;
    
    //Count the number of points
    if (!Count_points_in_cluster(cluster, structure, n_points)) {
        hout<<"Error in Export_cnts_in_cluster when calling Count_points_in_cluster"<<endl;
        return 0;
    }
    
    //Check that there are points
    if (!n_points) {
        hout<<"No points to export. NO visualization file was exported for CNTs."<<endl;
        return 1;
    }
    
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_cnts_in_cluster when calling Add_header"<<endl;
        return 0;
    }
    
    //Add the line with the number of points
    otec<<"POINTS "<<n_points<<" float" <<endl;
    
    //Add all the points
    if (!Add_points_in_cluster(points, structure, cluster, otec)) {
        hout<<"Error in Export_cnts_in_cluster when calling Add_points_in_cluster"<<endl;
        return 0;
    }
    
    //Add the line indicating the number of lines+1, and the number of points in those lines
    otec<<"LINES "<<cluster.size()+1<<' '<<n_points<<endl;
    
    //Add the offsets:
    //The number of points used after adding each line, starting with a zero
    if (!Add_offsets_from_cluster(structure, cluster, otec)) {
        hout<<"Error in Export_cnts_in_cluster when calling Add_offsets_from_cluster"<<endl;
        return 0;
    }
    
    //Add connectivity
    if (!Add_connectivity_from_cluster(structure, cluster, otec)) {
        hout<<"Error in Export_cnts_in_cluster when calling Add_connectivity_from_cluster"<<endl;
        return 0;
    }
    
    //Close the file
    otec.close();
    
    return 1;
}
int VTK_Export::Count_points_in_cluster(const vector<int> &cluster, const vector<vector<long int> > &structure, long int &n_points)const
{
    //Iterate over all lines (CNTs in the cluster)
    for (int i = 0; i < (int)cluster.size(); i++) {
        
        //Get the CNT number
        int CNTi = cluster[i];
        
        //Add the number of points in the CNT
        n_points = n_points + (long int)structure[CNTi].size();
    }
    
    return 1;
}
int VTK_Export::Add_points_in_cluster(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const vector<int> &cluster, ofstream &otec)const
{
    //Add the point coordinates, separated by spaces and in groups of four points
    //(12 coordinates) per line
    
    //Variable to count the points
    long int n_points = 1;
    
    //Iterate over all CNTs in the cluster
    for (int i = 0; i < (int)cluster.size(); i++) {
        
        //Get current CNT number
        int CNTi = cluster[i];
        
        //Iterate over all points in CNTi
        for (int j = 0; j < (int)structure[CNTi].size(); j++) {
            
            //Check if a new line needs to be started
            if (!(n_points%4)) {
                otec<<endl;
            }
            
            //Get point number
            long int P = structure[CNTi][j];
            
            //Append point i to file
            otec<<points[P].x<<' '<<points[P].y<<' '<<points[P].z<<' ';
            
            //Increase the count of points (for adding a new line only)
            n_points++;
        }
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
int VTK_Export::Add_offsets_from_cluster(const vector<vector<long int> > &structure, const vector<int> &cluster, ofstream &otec)const
{
    //Add the line with the OFFSETS command
    otec<<"OFFSETS vtktypeint64"<<endl;
    
    //Output a zero, which is required
    otec<<"0 ";
    
    //Accumulate the number of points
    long int acc_points = 0;
    
    //Count the number of offsets
    int n_offsets = 1;
    
    //Iterate over all CNTs in the cluster
    for (int i = 0; i < (int)cluster.size(); i++) {
        
        //Check if a new line needs to be started
        if (!(n_offsets%20)) {
            otec<<endl;
        }
        
        //Get current CNT number
        int CNTi = cluster[i];
        
        //Accumulate the number of points
        acc_points = acc_points + (long int)structure[CNTi].size();
        
        //Output the accumulated number of points
        otec<<acc_points<<' ';
        
        //Each element in cluster is a line (CNT), so increase the number of lines
        n_offsets++;
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
int VTK_Export::Add_connectivity_from_cluster(const vector<vector<long int> > &structure, const vector<int> &cluster, ofstream &otec)const
{
    //Add the line with the CONNECTIVITY command
    otec<<"CONNECTIVITY vtktypeint64"<<endl;
    
    //Variable to count the number of points
    long int n_points = 0;
    
    //Add the connectivity, one connectivity per line
    //Iterate over all CNTs in the cluster
    for (int i = 0; i < (int)cluster.size(); i++) {
        
        //Get current CNT number
        int CNTi = cluster[i];
        
        //Iterate over all points in CNTi
        for (int j = 0; j < (int)structure[CNTi].size(); j++) {
            
            //Add the consecutive number of point structure[CNTi][j] in a line
            otec<<n_points<<' ';
            
            //Increase the number of points
            n_points++;
        }
        
        //For every line (CNT) start a new line in the file
        otec<<endl;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Function to export a single CNT
int VTK_Export::Export_single_cnt(const vector<Point_3D> &points, const string &filename)const
{
    //Check that there are at least two points in the CNT
    //(a single point does not make a CNT)
    if (points.size() <= 1) {
        hout<<"No point in the CNT to export. NO visualization file was exported. Input filename: "<<filename<<endl;
        return 1;
    }
    
    //Count the number of points
    int n_points = (int)points.size();
    
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_single_cnt when calling Add_header"<<endl;
        return 0;
    }
    
    //Add the line with the number of points
    otec<<"POINTS "<<n_points<<" float" <<endl;
    
    //Add all the points
    if (!Add_points_from_vector(points, otec)) {
        hout<<"Error in Export_single_cnt when calling Add_points_from_vector"<<endl;
        return 0;
    }
    
    //Add the line indicating number of lines+1 (i.e., 2), and the number of points in that line
    otec<<"LINES 2 "<<n_points<<endl;
    
    //Add the offsets
    //Add the line with the OFFSETS command
    otec<<"OFFSETS vtktypeint64"<<endl;
    
    //First output a zero, which is required then output the number of points in the line (CNT)
    otec<<"0 "<<n_points<<endl;
    
    //Add connectivity
    //Add the line with the CONNECTIVITY command
    otec<<"CONNECTIVITY vtktypeint64"<<endl;
    
    //Iterate over all points in the CNT
    otec<<"0 ";
    for (size_t j = 1; j < points.size(); j++) {
        
        //Add the consecutive number of point j within the CNT
        otec<<j<<' ';
        
        //Add a new line every 50 points
        if ( !(j%50) ) {
            otec<<endl;
        }
    }
    //Add one line at the end
    otec<<endl;
    
    //Close the file
    otec.close();
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Function to export CNTs from the vector of points and the structure
int VTK_Export::Export_from_cnt_structure(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const string &filename)const
{
    //Check that there are CNTs to export
    if (points.empty() || structure.empty()) {
        hout<<"No CNTs to export. NO visualization file was exported for CNTs."<<endl;
        return 1;
    }
    
    //Variable to store the number of points
    long int n_points = 0;
    
    //Variable to store the number of lines (CNTs)
    int n_cnts = (int)structure.size();
    
    //Count the number of points
    if (!Count_points_in_structure(structure, n_points)) {
        hout<<"Error in Export_from_cnt_structure when calling Count_points_in_structure"<<endl;
        return 0;
    }
    
    //Check that the points vector has points
    if (!n_points) {
        hout<<"No points to export. NO visualization file was exported for CNTs."<<endl;
        return 1;
    }
    
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_from_cnt_structure when calling Add_header"<<endl;
        return 0;
    }
    
    //Add the line with the number of points
    otec<<"POINTS "<<n_points<<" float" <<endl;
    
    //Add all the points
    if (!Add_points_from_structure(points, structure, otec)) {
        hout<<"Error in Export_from_cnt_structure when calling Add_points_from_structure"<<endl;
        return 0;
    }
    
    //Add the line indicating the number of lines+1, and the number of points in those lines
    otec<<"LINES "<<n_cnts+1<<' '<<n_points<<endl;
    
    //Add the offsets:
    //The number of points used after adding each line, starting with a zero
    if (!Add_offsets_from_structure(structure, otec)) {
        hout<<"Error in Export_from_cnt_structure when calling Add_offsets_from_structure"<<endl;
        return 0;
    }
    
    //Add connectivity
    if (!Add_connectivity_from_structure(structure, otec)) {
        hout<<"Error in Export_from_cnt_structure when calling Add_connectivity_from_structure"<<endl;
        return 0;
    }
    
    //Close the file
    otec.close();
    
    return 1;
}
int VTK_Export::Count_points_in_structure(const vector<vector<long int> > &structure, long int &n_points)const
{
    //Iterate over all lines (CNTs in the structure)
    for (int i = 0; i < (int)structure.size(); i++) {
        
        //Add the number of points in the CNT
        n_points = n_points + (long int)structure[i].size();
    }
    
    return 1;
}
//This function adds the coordinates of a vector of points into a file as indicated
//by the structure vector
//Four points are printed per line
int VTK_Export::Add_points_from_structure(const vector<Point_3D> &points, const vector<vector<long int> > &structure, ofstream &otec)const
{
    //Add the point coordinates, separated by spaces and in groups of four points
    //(12 coordinates) per line
    
    //Variable to count the points
    long int n_points = 1;
    
    //Iterate over all the pairs of indices
    for (int i = 0; i < (int)structure.size(); i++) {
        for (int j = 0; j < (int)structure[i].size(); j++) {
            
            //Check if a new line needs to be started
            if (!(n_points%4)) {
                otec<<endl;
            }
            
            //Get current point of the line (CNT)
            long int P1 = structure[i][j];
            
            //Append point i to file
            otec<<points[P1].x<<' '<<points[P1].y<<' '<<points[P1].z<<' ';
            
            //Increase the count of points (for adding a new line only)
            n_points++;
        }
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
int VTK_Export::Add_offsets_from_structure(const vector<vector<long int> > &structure, ofstream &otec)const
{
    //Add the line with the OFFSETS command
    otec<<"OFFSETS vtktypeint64"<<endl;
    
    //Output a zero, which is required
    otec<<"0 ";
    
    //Accumulate the number of points
    long int acc_points = 0;
    
    //Count the number of offsets
    int n_offsets = 0;
    
    //Iterate over all the pairs of indices
    for (int i = 0; i < (int)structure.size(); i++) {
        
        //Add the difference of indices plus one
        acc_points = acc_points + (long int)structure[i].size();
        
        //Output the accumulated number of points
        otec<<acc_points<<' ';
        
        //Increase the number of lines
        n_offsets++;
        
        //Check if a new line needs to be started
        if (!(n_offsets%20)) {
            otec<<endl;
        }
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
int VTK_Export::Add_connectivity_from_structure(const vector<vector<long int> > &structure, ofstream &otec)const
{
    //Add the line with the CONNECTIVITY command
    otec<<"CONNECTIVITY vtktypeint64"<<endl;
    
    //Variable to count the number of points
    long int n_points = 0;
    
    //Add the connectivity, one connectivity per line
    for (int i = 0; i < (int)structure.size(); i++) {
        for (int j = 0; j < (int)structure[i].size(); j++) {
            
            //Add the consecutive number of point structure[i][j] in a line
            otec<<n_points<<' ';
            
            //Increase the number of points
            n_points++;
            
            //Check if a new line needs to be started
            if (!(n_points%50)) {
                otec<<endl;
            }
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//GNPs
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
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
    if (!Add_ofsets_for_gnps((int)gnps.size(), otec)) {
        hout<<"Error in Export_gnps when calling Add_ofsets_for_gnps"<<endl;
        return 0;
    }
    
    //Add the connectivity
    if (!Add_connectivity_for_gnps((long int)gnps.size(), otec)) {
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
    otec.precision(15);
    
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
int VTK_Export::Add_ofsets_for_gnps(const int &n_gnps, ofstream &otec)const
{
    //Add the line with the OFFSETS command
    otec<<"OFFSETS vtktypeint64"<<endl;
    
    //Output a zero, which is required
    otec<<"0 ";
    
    //Variable to count the number of vertices
    int n_vertices = 0;
    
    //Write the offsets, one line per GNP
    for (size_t i = 0; i < n_gnps; i++) {
        
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
int VTK_Export::Add_connectivity_for_gnps(const long int &n_gnps, ofstream &otec, const long int &cnt_points_offset)const
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
    
    for (long int i = 0; i < n_gnps; i++) {
        
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
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int VTK_Export::Export_gnps_in_cluster(const vector<GNP> &gnps, const vector<int> &cluster, const string &filename)
{
    //Check that the points vector has points
    if (cluster.empty()) {
        hout<<"Vector of GNPs is empty. NO visualization file was exported."<<endl;
        return 1;
    }
    
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_gnps_in_cluster when calling Add_header"<<endl;
        return 0;
    }
    
    //Add the line with the number of points, which is the number of GNPs multiplied by 8
    otec<<"POINTS "<<cluster.size()*8<<" float" <<endl;
    
    //Add all the points, i.e., all the vertices of the GNPs
    if (!Add_all_gnp_vertices_from_cluster(gnps, cluster, otec)) {
        hout<<"Error in Export_gnps_in_cluster when calling Add_all_gnp_vertices_from_cluster"<<endl;
        return 0;
    }
    
    //Add the line with the polygons command, with the number of faces+1 and
    //the number of points in those faces
    otec<<"POLYGONS "<<cluster.size()*6 + 1<<' '<<cluster.size()*24<<endl;
    
    //Add the offsets
    if (!Add_ofsets_for_gnps((int)cluster.size(), otec)) {
        hout<<"Error in Export_gnps_in_cluster when calling Add_ofsets_for_gnps"<<endl;
        return 0;
    }
    
    //Add the connectivity
    if (!Add_connectivity_for_gnps((long int)cluster.size(), otec)) {
        hout<<"Error in Export_gnps_in_cluster when calling Add_connectivity_for_gnps"<<endl;
        return 0;
    }
    
    //Close the file
    otec.close();
    
    return 1;
}
int VTK_Export::Add_all_gnp_vertices_from_cluster(const vector<GNP> &gnps, const vector<int> &cluster, ofstream &otec)const
{
    //Iterate over the GNPs in the cluster
    for (size_t i = 0; i < cluster.size(); i++) {
        
        //Get the GNP number
        int GNPi = cluster[i];
        
        //Add the vertices
        if (!Add_points_from_array(gnps[GNPi].vertices, 8, otec)) {
            hout<<"Error in Add_all_gnp_vertices_from_cluster when calling Add_points_from_array"<<endl;
            return 0;
        }
        
        //Add a new line
        otec<<endl;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//This function exports all GNPs in a vector in single files, i.e., one GNP per visualization file
int VTK_Export::Export_gnps_single_files(const vector<GNP>& gnps, const string& base_filename)const
{
    //Check that the points vector has points
    if (gnps.empty()) {
        hout << "Vector of GNPs is empty. NO visualization file was exported." << endl;
        return 1;
    }

    //Iterate over the GNPs in the vector
    for (size_t i = 0; i < gnps.size(); i++)
    {
        //Generate file name
        string filename = base_filename + "_" + to_string(i) + ".vtk";
         
        //Export GNP i
        if (!Export_single_gnp(gnps[i], filename))
        {
            hout << "Error in Export_gnps_single_files when calling Export_single_gnp for iteration " << i << " of " << gnps.size() << endl;
            return 0;
        }
    }

    return 1;
}
//This function exports a single GNP into a visualization file
int VTK_Export::Export_single_gnp(const GNP& gnp_i, const string& filename)const
{
    //Open the file
    ofstream otec(filename.c_str());

    //Add header
    if (!Add_header(otec)) {
        hout << "Error in Export_single_gnp when calling Add_header" << endl;
        return 0;
    }

    //Add the line with the number of points, which is the case of a single GNP is 8
    otec << "POINTS " <<  8 << " float" << endl;

    //Add the vertices
    if (!Add_points_from_array(gnp_i.vertices, 8, otec)) {
        hout << "Error in Export_single_gnp when calling Add_points_from_array" << endl;
        return 0;
    }

    //Add a new line
    otec << endl;

    //Add the line with the polygons command, with the number of faces+1 and
    //the number of points in those faces
    otec << "POLYGONS " << 7 << ' ' << 24 << endl;

    //Add the offsets
    if (!Add_ofsets_for_gnps(1, otec)) {
        hout << "Error in Export_single_gnp when calling Add_ofsets_for_gnps" << endl;
        return 0;
    }

    //Add the connectivity
    if (!Add_connectivity_for_gnps(1, otec)) {
        hout << "Error in Export_single_gnp when calling Add_connectivity_for_gnps" << endl;
        return 0;
    }

    //Close the file
    otec.close();


    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Cuboid
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int VTK_Export::Export_cuboid(const cuboid &cub, const string &filename)
{
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_gnps when calling Add_header"<<endl;
        return 0;
    }
    
    //Add the line with the number of points in the cuboid
    otec<<"POINTS 8 float" <<endl;
    
    //Add all the points, i.e., all the vertices of the cuboid
    otec<<cub.poi_min.x<<' '<<cub.poi_min.y<<' '<<cub.poi_min.z<<' ';
    otec<<cub.poi_min.x+cub.len_x<<' '<<cub.poi_min.y<<' '<<cub.poi_min.z<<' ';
    otec<<cub.poi_min.x+cub.len_x<<' '<<cub.poi_min.y+cub.wid_y<<' '<<cub.poi_min.z<<' ';
    otec<<cub.poi_min.x<<' '<<cub.poi_min.y+cub.wid_y<<' '<<cub.poi_min.z<<' ';
    otec<<endl;
    otec<<cub.poi_min.x<<' '<<cub.poi_min.y<<' '<<cub.poi_min.z+cub.hei_z<<' ';
    otec<<cub.poi_min.x+cub.len_x<<' '<<cub.poi_min.y<<' '<<cub.poi_min.z+cub.hei_z<<' ';
    otec<<cub.poi_min.x+cub.len_x<<' '<<cub.poi_min.y+cub.wid_y<<' '<<cub.poi_min.z+cub.hei_z<<' ';
    otec<<cub.poi_min.x<<' '<<cub.poi_min.y+cub.wid_y<<' '<<cub.poi_min.z+cub.hei_z<<endl;
    
    //Add the line with the polygons command, with the number of faces+1 and
    //the number of points in those faces
    otec<<"POLYGONS 7 24"<<endl;
    
    //Add the line with the OFFSETS command
    otec<<"OFFSETS vtktypeint64"<<endl;
    
    //Add the offsets
    otec<<"0 4 8 12 16 20 24"<<endl;
    
    //Add the line with the CONNECTIVITY command
    otec<<"CONNECTIVITY vtktypeint64"<<endl;
    
    //Add the connectivity
    otec<<"0 1 2 3"<<endl;
    otec<<"4 5 6 7"<<endl;
    otec<<"3 0 4 7"<<endl;
    otec<<"0 1 5 4"<<endl;
    otec<<"1 2 6 5"<<endl;
    otec<<"2 3 7 6"<<endl;
    
    //Close the file
    otec.close();
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Mixed
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
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
    if (!Add_ofsets_for_gnps((int)gnps.size(), otec)) {
        hout<<"Error in Export_hybrid_material when calling Add_ofsets_for_gnps"<<endl;
        return 0;
    }
    
    //Add the connectivity
    //the first GNP point number is points.size(), i.e., the offset is points.size()
    if (!Add_connectivity_for_gnps((long int)gnps.size(), otec, (long int)points.size())) {
        hout<<"Error in Export_hybrid_material when calling Add_connectivity_for_gnps"<<endl;
        return 0;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Triangulations
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int VTK_Export::Export_triangulation(const vector<Point_3D> &points, const vector<EdgeL> &triangulation, const string &filename)const
{
    //Check that the traingulation has edges
    if (triangulation.empty()) {
        hout<<"Vector of triangulation edges is empty. NO visualization file was exported."<<endl;
        return 1;
    }
    
    //Open the file
    ofstream otec(filename.c_str());
    
    //Add header
    if (!Add_header(otec)) {
        hout<<"Error in Export_triangulation when calling Add_header"<<endl;
        return 0;
    }
    
    //Add the line with the number of points, which is twice the number of edges
    otec<<"POINTS "<<2*(triangulation.size())<<" float" <<endl;
    
    //Add all the points
    if (!Add_points_from_triangulation_edges(points, triangulation, otec)) {
        hout<<"Error in Export_triangulation when calling Add_points_from_triangulation_edges"<<endl;
        return 0;
    }
    
    //Add the line indicating the number of lines+1, and the number of points in those lines
    otec<<"LINES "<<triangulation.size()+1<<' '<<2*(triangulation.size())<<endl;
    
    //Add the offsets:
    //The number of points used after adding each line, starting with a zero
    if (!Add_offsets_for_trinagulation(triangulation, otec)) {
        hout<<"Error in Export_triangulation when calling Add_offsets_for_trinagulation"<<endl;
        return 0;
    }
    
    //Add connectivity
    if (!Add_connectivity_for_trinagulation(triangulation, otec)) {
        hout<<"Error in Export_triangulation when calling Add_connectivity_for_trinagulation"<<endl;
        return 0;
    }
    
    //Close the file
    otec.close();
    
    return 1;
}
int VTK_Export::Add_points_from_triangulation_edges(const vector<Point_3D> &points, const vector<EdgeL> &triangulation, ofstream &otec)const
{
    
    //Add the point coordinates, separated by spaces and in groups of four points
    //(12 coordinates) per line
    
    //Add the points from the first vertex
    long int v = triangulation[0].v1;
    otec<<points[v].x<<' '<<points[v].y<<' '<<points[v].z<<' ';
    v = triangulation[0].v2;
    otec<<points[v].x<<' '<<points[v].y<<' '<<points[v].z<<' ';
    
    //Add the remaninig points
    for (size_t i = 1; i < triangulation.size(); i++) {
        
        //Check if a new line needs to be started
        if (!(i%4)) {
            otec<<endl;
        }
        
        //Add points of edge i
        v = triangulation[i].v1;
        otec<<points[v].x<<' '<<points[v].y<<' '<<points[v].z<<' ';
        v = triangulation[i].v2;
        otec<<points[v].x<<' '<<points[v].y<<' '<<points[v].z<<' ';
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
int VTK_Export::Add_offsets_for_trinagulation(const vector<EdgeL> &triangulation, ofstream &otec)const
{
    //Add the line with the OFFSETS command
    otec<<"OFFSETS vtktypeint64"<<endl;
    
    //Output a zero, which is required
    otec<<"0 ";
    
    //Variable to store the accumulated number of points,
    //initialize with 2 which is the number of points in the first vertex
    long int acc_points = 2;
    
    //Output the accumulated number of points
    otec<<acc_points<<' ';
    
    //Output the accumulated number of points for the remaining edges
    //Print 20 per line
    for (size_t i = 1; i < triangulation.size(); i++) {
        
        //Check if a new line needs to be started
        if (!(i%20)) {
            otec<<endl;
        }
        
        //Add the number of points from edge i to the accumulated number of points
        //i.e., 2 points in each edge
        acc_points += 2;
        
        //Output the accumulated number of points
        otec<<acc_points<<' ';
    }
    
    //Start a new line
    otec<<endl;
    
    return 1;
}
int VTK_Export::Add_connectivity_for_trinagulation(const vector<EdgeL> &triangulation, ofstream &otec)const
{
    //Add the line with the CONNECTIVITY command
    otec<<"CONNECTIVITY vtktypeint64"<<endl;
    
    //Variable to count the number of points
    long int n_points = 0;
    
    //Add the connectivity, one connectivity per line
    for (size_t i = 0; i < triangulation.size(); i++) {
        
        //Add the consecutive number of points in the vertex
        //Not the vertex number, but the consecutive number in which they were printed
        //in th evtk file
        otec<<n_points<<' '<<n_points+1<<endl;
        
        //Increase the count of points
        n_points = n_points + 2;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Triangulations
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int VTK_Export::Export_triangles(const vector<Point_3D>& vertices, const vector<TrFaceL>& triangles, const string& filename)const
{
    //Check that the traingulation has edges
    if (triangles.empty()) {
        hout << "Vector of triangles edges is empty. NO visualization file was exported." << endl;
        return 1;
    }

    //Open the file
    ofstream otec(filename.c_str());

    //Add header
    if (!Add_header(otec)) {
        hout << "Error in Export_triangles when calling Add_header" << endl;
        return 0;
    }

    //Add the line with the number of points, which is twice the number of vertices
    otec << "POINTS " << vertices.size() << " float" << endl;

    //Add all the points
    if (!Add_points_from_vector(vertices, otec)) {
        hout << "Error in Export_triangles when calling Add_points_from_vector" << endl;
        return 0;
    }

    //Get the number of triangles
    int n_tri = (int)triangles.size();

    //Add the line with the polygons command, with the number of faces+1 and
    //the number of points in those faces
    otec << "POLYGONS " << (n_tri + 1) << ' ' << 3*n_tri << endl;

    //Add the offsets
    if (!Add_ofsets_for_triangles(n_tri, otec)) {
        hout << "Error in Export_triangles when calling Add_ofsets_for_triangles" << endl;
        return 0;
    }

    //Add the connectivity
    if (!Add_connectivity_for_triangles(triangles, otec)) {
        hout << "Error in Export_single_gnp when calling Add_connectivity_for_triangles" << endl;
        return 0;
    }

    //Close the file
    otec.close();

    return 1;
}
//
int VTK_Export::Add_ofsets_for_triangles(const int& n_tri, ofstream& otec)const
{
    //Add the line with the OFFSETS command
    otec << "OFFSETS vtktypeint64" << endl;

    //Output a zero, which is required
    otec << "0 ";

    //Variable to store the accumulated number of points
    //Initialized with 3, which is the number of points (vertices) in a triangle
    int acc_pts = 3;

    //Iterate over the number of triangles
    for (int i = 1; i <= n_tri; i++)
    {
        //Output the number of accumulated points
        otec << acc_pts << " ";

        //Increase the accumulated number of points
        acc_pts = acc_pts + 3;

        //Check if a new line is added
        if (!(i%20))
        {
            otec << endl;
        }
    }

    //Add a new line
    otec << endl;

    return 1;
}
int VTK_Export::Add_connectivity_for_triangles(const vector<TrFaceL>& triangles, ofstream& otec)const
{
    //Add the line with the CONNECTIVITY command
    otec << "CONNECTIVITY vtktypeint64" << endl;

    //Iterate over the triangles
    for (size_t i = 0; i < triangles.size(); i++)
    {
        //Export vertices of triangle i
        otec << triangles[i].v1 << " " << triangles[i].v2 << " " << triangles[i].v3 << " ";

        //Check if a new line is added
        if (!(i % 10))
        {
            otec << endl;
        }
    }

    //Add a new line
    otec << endl;

    return 1;
}
int VTK_Export::Export_supertetrahedron(const vector<Point_3D>& vertices, const vector<TrFaceL>& triangles, const string& filename)const
{
    //Check that the traingulation has edges
    if (triangles.empty()) {
        hout << "Vector of supertetrahedron triangles edges is empty. NO visualization file was exported." << endl;
        return 1;
    }

    //Generate a new vector of triangles the same size as triangles
    vector<TrFaceL> new_triangles(triangles.size());

    //Transform all vertices of the supertetrahedron into valid indices
    for (size_t i = 0; i < triangles.size(); i++)
    {
        new_triangles[i].v1 = -triangles[i].v1 - 1;
        new_triangles[i].v2 = -triangles[i].v2 - 1;
        new_triangles[i].v3 = -triangles[i].v3 - 1;
    }

    //Export triangles using the new vector of triangles
    if (!Export_triangles(vertices, new_triangles, filename))
    {
        hout << "Error in Export_supertetrahedron when calling Export_triangles" << endl;
        return 0;
    }

    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Points
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//This function exports an array of points
int VTK_Export::Export_point_array(const Point_3D points[], const int& size, const string& filename)const
{
    //Open the file
    ofstream otec(filename.c_str());

    //Add header
    if (!Add_header(otec)) 
    {
        hout << "Error in Export_point_array when calling Add_header" << endl;
        return 0;
    }

    //Add point coordinates
    if (!Add_point_coordinates_from_array(points, size, otec))
    {
        hout << "Error in Export_point_array when calling Add_point_coordinates_from_array" << endl;
        return 0;
    }

    //Add the lines corresponding to vertices, offsets and connectivity
    if (!Add_vertices_offsets_connectivity_for_n_points(size, otec))
    {
        hout << "Error in Export_point_array when calling Add_vertices_offsets_connectivity_for_n_points" << endl;
    }

    //Close the file
    otec.close();

    return 1;
}
//This function adds the number of points and the point coordinates when exporting a point array
int VTK_Export::Add_point_coordinates_from_array(const Point_3D points[], const int& size, ofstream& otec)const
{
    //Add the number of points
    otec << "POINTS " << size << " float";

    //Add the point coordinates
    for (int i = 0; i < size; i++)
    {
        //Add a line break every four points
        if (!(i % 4))
            otec << endl;
        
        //Add point coordinates separated by a space
        otec << points[i].x << " " << points[i].y << " " << points[i].z << " ";
    }

    //Add a final line break
    otec << endl;

    return 1;
}
//
int VTK_Export::Add_vertices_offsets_connectivity_for_n_points(const int& n, ofstream& otec)const
{
    //Add vertices line
    //The first number is the number of points+1 and the second number is the number of points
    otec << "VERTICES" << n + 1 << " " << n << endl;

    //Add offsets line
    otec << "OFFSETS vtktypeint64";

    //Add offsets (sequence up to number of points)
    if (!Add_consecutive_numbers(n, otec))
    {
        hout << "Error in Add_vertices_offsets_connectivity_for_n_points when calling Add_consecutive_numbers (offsets)" << endl;
        return 0;
    }

    //Add connectivity line
    otec << "CONNECTIVITY vtktypeint64";

    //Add connectivity (sequence up to number of points - 1)
    if (!Add_consecutive_numbers(n - 1, otec))
    {
        hout << "Error in Add_vertices_offsets_connectivity_for_n_points when calling Add_consecutive_numbers (connectivity)" << endl;
        return 0;
    }

    return 1;
}
//This function adds consecutive numbers from 0 to n (including n)
int VTK_Export::Add_consecutive_numbers(const int& n, ofstream& otec)const
{
    //Iterate from 0 to n (including n)
    for (int i = 0; i <= n; i++)
    {
        //Add a line break every 20 points
        if (!(i % 20))
            otec << endl;

        //Add offset
        otec << i << " ";
    }

    //Add a final line break
    otec << endl;

    return 1;
}
//This function exports a vector of points as points in a VTK file
int VTK_Export::Export_point_vector(const vector<Point_3D>& points, const string& filename)const
{
    //Check that the vector is not empty
    if (points.empty()) {
        hout << "Vector of points is empty. NO visualization file was exported." << endl;
        return 1;
    }

    //Open the file
    ofstream otec(filename.c_str());

    //Add header
    if (!Add_header(otec))
    {
        hout << "Error in Export_point_vector when calling Add_header" << endl;
        return 0;
    }

    //Add point coordinates
    if (!Add_point_coordinates_from_vector(points, otec))
    {
        hout << "Error in Export_point_vector when calling Add_point_coordinates_from_vector" << endl;
        return 0;
    }

    //Get the number of points
    int n_p = (int)points.size();

    //Add the lines corresponding to vertices, offsets and connectivity
    if (!Add_vertices_offsets_connectivity_for_n_points(n_p, otec))
    {
        hout << "Error in Export_point_vector when calling Add_vertices_offsets_connectivity_for_n_points" << endl;
    }

    //Close the file
    otec.close();

    return 1;
}
//
int VTK_Export::Add_point_coordinates_from_vector(const vector<Point_3D>& points, ofstream& otec)const
{
    //Get the number of points
    int n_p = (int)points.size();

    //Add the number of points
    otec << "POINTS " << n_p << " float";

    //Add the point coordinates
    for (int i = 0; i < n_p; i++)
    {
        //Add a line break every four points
        if (!(i % 4))
            otec << endl;

        //Add point coordinates separated by a space
        otec << points[i].x << " " << points[i].y << " " << points[i].z << " ";
    }

    //Add a final line break
    otec << endl;

    return 1;
}
//This function exports the points in an array but only those indicated in a vector of vertices
int VTK_Export::Export_selected_points_in_array(const vector<int>& vertices, const Point_3D points[], const string& filename)const
{
    //Check that the vector of vertices is not empty
    if (vertices.empty()) {
        hout << "Vector of vertices to export is empty. NO visualization file was exported." << endl;
        return 1;
    }

    //Open the file
    ofstream otec(filename.c_str());

    //Add header
    if (!Add_header(otec))
    {
        hout << "Error in Export_selected_points_in_array when calling Add_header" << endl;
        return 0;
    }

    //Add point coordinates
    if (!Add_point_coordinates_from_vertex_vector(vertices, points, otec))
    {
        hout << "Error in Export_selected_points_in_array when calling Add_point_coordinates_from_vertex_vector" << endl;
        return 0;
    }
    //Get the number of points
    int n_p = (int)vertices.size();

    //Add the lines corresponding to vertices, offsets and connectivity
    if (!Add_vertices_offsets_connectivity_for_n_points(n_p, otec))
    {
        hout << "Error in Export_point_vector when calling Add_vertices_offsets_connectivity_for_n_points" << endl;
    }

    //Close the file
    otec.close();

    return 1;
}

int VTK_Export::Add_point_coordinates_from_vertex_vector(const vector<int>& vertices, const Point_3D points[], ofstream& otec)const
{
    //Get the number of points
    int n_p = (int)vertices.size();

    //Add the number of points
    otec << "POINTS " << n_p << " float";

    //Add the point coordinates
    for (int i = 0; i < n_p; i++)
    {
        //Add a line break every four points
        if (!(i % 4))
            otec << endl;

        //Get the vertex number
        int v = vertices[i];

        //Add point coordinates separated by a space
        otec << points[v].x << " " << points[v].y << " " << points[v].z << " ";
    }

    //Add a final line break
    otec << endl;

    return 1;
}
