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

//Print the four vertices of a GNP needed to generated them in Abaqus
void Printer::Print_4_vertices_gnps(const vector<GNP> &gnps, const string &filename)
{
    //Open file
    ofstream otec(filename.c_str());
    otec.precision(10);
    
    //Iterate over all GNPs
    for (size_t i = 0; i < gnps.size(); i++) {
        
        //Ouput the coordinates of vertices 0 to 2 and 4
        otec<<gnps[i].vertices[0].x<<' '<<gnps[i].vertices[0].y<<' '<<gnps[i].vertices[0].z<<' '<<endl;
        otec<<gnps[i].vertices[1].x<<' '<<gnps[i].vertices[1].y<<' '<<gnps[i].vertices[1].z<<' '<<endl;
        otec<<gnps[i].vertices[2].x<<' '<<gnps[i].vertices[2].y<<' '<<gnps[i].vertices[2].z<<' '<<endl;
        otec<<gnps[i].vertices[4].x<<' '<<gnps[i].vertices[4].y<<' '<<gnps[i].vertices[4].z<<' '<<endl;
    }
    
    //Close file
    otec.close();
}
