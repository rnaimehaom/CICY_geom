//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Functions to print (output into a file) different types of data structures
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Printer.h"


//Print a vector of 3D points
void Printer::Print_1d_vec(const vector<Point_3D> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    otec.precision(15);
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)list.size(); i++) {
        otec << list[i].x << "\t" << list[i].y << "\t" << list[i].z << "\t" << list[i].flag << "\n";
    }
    otec.close();
}

//Print a vector of chars with the specified filename
void Printer::Print_1d_vec(const vector<char> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of integers with the specified filename
void Printer::Print_1d_vec(const vector<int> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    //otec << "Positions \n";
    for (long int i=0; i < (long int)list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of doubles with the specified filename
void Printer::Print_1d_vec(const vector<double> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    otec.precision(15);
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of doubles with the specified filename
void Printer::Append_1d_vec(const vector<double> &list, const string &filename)
{
    ofstream otec(filename.c_str(), std::ios_base::app);
    otec.precision(15);
    hout << "Appending to file: " << filename << "\n";
    for (long int i=0; i < (long int)list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\t";
    }
    otec << "\n";
    otec.close();
}

//Print a vector of long integers with the specified filename
void Printer::Print_1d_vec(const vector<long int> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of vectors of integers with the specified filename
void Printer::Print_2d_vec(const vector<vector<int> > &num_mat, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)num_mat.size(); i++) {
        for (long int j = 0; j < (long int)num_mat[i].size(); j++) {
            otec << num_mat[i][j] << '\t' ;
        }
        otec << '\n' ;
    }
    otec.close();
}

//Print a vector of vectors of integers with the specified filename
void Printer::Print_2d_vec(const vector<vector<long int> > &num_mat, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)num_mat.size(); i++) {
        for (long int j = 0; j < (long int)num_mat[i].size(); j++) {
            otec << num_mat[i][j] << '\t' ;
        }
        otec << '\n' ;
    }
    otec.close();
}

//Print a vector of vectors of doubles with the specified filename
void Printer::Print_2d_vec(const vector<vector<double> > &num_mat, const string &filename)
{
    ofstream otec(filename.c_str());
    otec.precision(15);
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)num_mat.size(); i++) {
        for (long int j = 0; j < (long int)num_mat[i].size(); j++) {
            otec << num_mat[i][j] << '\t' ;
        }
        otec << '\n' ;
    }
    otec.close();
}

//Print the GNP data needed to generate them in Abaqus
void Printer::Print_gnp_data(const vector<GNP> &gnps, const int &prec, const string &filename)
{
    //Open file
    ofstream otec(filename.c_str());
    
    //Calculate 2PI so the operation is done only once
    double two_PI = 2*PI;
    
    //Iterate over all GNPs
    for (size_t i = 0; i < gnps.size(); i++) {
        
        //Output the geometry of the GNP separated by commas
        otec<<gnps[i].l<<", "<<gnps[i].l<<", "<<gnps[i].t<<", ";
        
        //Calculate the rotation angles
        //Angle phi corresponds to the rotation angle around z
        double phi = Recover_angle(gnps[i].rotation.element[1][1], abs(gnps[i].rotation.element[0][1]), two_PI);
        //Angle theta corresponds to the rotation angle around y
        double theta = Recover_angle(gnps[i].rotation.element[2][2], abs(gnps[i].rotation.element[2][0]), two_PI);
        
        //Output the rotation angles
        otec<<theta<<", "<<phi<<", ";
        
        //Ouput the coordinates of the centroid and add a line break
        otec<<gnps[i].center.str(prec)<<endl;
    }
    
    //Close file
    otec.close();
}

//This function is used to recover the two angles used in the rotation matrix
double Printer::Recover_angle(const double &cos_alpha, const double &sin_alpha, const double &two_PI)
{
    
    //Check the signs of the trigonometric functions
    if (cos_alpha > Zero) {
        if (sin_alpha > Zero) {
            
            //No adjustment required
            return acos(cos_alpha);
        }
        else {
            
            //Recover angle from the sine
            double phi = asin(sin_alpha) + two_PI;
            return phi;
        }
    }
    else {
        if (sin_alpha > Zero) {
            
            //Angle can be ontained directly fom cosine
            return acos(cos_alpha);
        }
        else {
            
            //Recover angle from the cosine
            double phi = two_PI - acos(cos_alpha);
            return phi;
        }
    }
}

//Print the four vertices of a GNP needed to generate them in Abaqus
//The precision (number of digits after the decimal point) is specified as an input
void Printer::Print_4_vertices_gnps(const vector<GNP> &gnps, const int &prec, const string &filename)
{
    //Open file
    ofstream otec(filename.c_str());
    
    //Iterate over all GNPs
    for (size_t i = 0; i < gnps.size(); i++) {
        
        //Ouput the coordinates of vertices 0 to 2 and 4
        otec<<gnps[i].vertices[0].str(prec)<<endl;
        otec<<gnps[i].vertices[1].str(prec)<<endl;
        otec<<gnps[i].vertices[2].str(prec)<<endl;
        otec<<gnps[i].vertices[4].str(prec)<<endl;
    }
    
    //Close file
    otec.close();
}

//This function prints the coordinates of all CNT points into a file
//Two files are exported, one with the coordinates and one with the number of CNTs and
//the number of points for each CNT
void Printer::Print_cnt_points_and_structure(const vector<vector<long int> > &structure, const vector<Point_3D> &points_cnt, const vector<double> &radii, const int &prec, const string &filename_points, const string &filename_struct)
{
    //Open file for CNT points
    ofstream otec_points(filename_points.c_str());
    //Open file for CNT structure
    ofstream otec_struct(filename_struct.c_str());
    
    //Output the number of CNTs and radii
    otec_struct<<structure.size()<<", "<<radii.size()<<endl;
    
    //Iterate over all points in the structure
    for (size_t i = 0; i < structure.size(); i++) {
        
        //Output the number of points in CNT i and its radius
        otec_struct<<structure[i].size()<<", "<<radii[i]<<endl;
        
        //Iterate over the points in CNT i
        for (size_t j = 0; j < structure[i].size(); j++) {
            
            //Get the point number
            long int Pj = structure[i][j];
            
            //Output the coordinates of point j in CNT i
            otec_points<<points_cnt[Pj].str(prec)<<endl;
        }
    }
    
    
    //Close files
    otec_points.close();
    otec_struct.close();
}
