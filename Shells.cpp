//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Generate shells (background vectors) to map each CNT into an observation window and facilitate their trimming
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Shells.h"

/*
 
 This class generates the vector:
 shells_cnt: is used to trim the CNTs. The CNTs are assigned to a shell sub-region (for them to be trimmed).
 The shells will be created as follows: The smallest observation window (cuboid) will be one shell sub-region and will be used for the last element in the vector.
 Then, the next shell-subregion will be the volume of the next observation window minus the volume of the first one (the smallest observation window).
 The following shells will have the same form: the volume of the observation window minus the volume of the previous one.
 The last shell region will be the boundary layer. This will be used for the first element of the vector
 
 */

int Shells::Generate_shells(const struct Geom_sample &sample, const vector<Point_3D> &points_cnt, const vector<GNP> &gnps, vector<vector<int> > &shells_cnt, vector<Shell> &shells_gnps)
{
    //Calculate the number of shells
    //sample.cut_num is the number of increments from the smallest to the largest observation widows
    //Then, there are sample.cut_num+1 observation windows
    //Then, there are sample.cut_num+2 shell sub-regions if we include the boundary layer.
    //Hence the shells_cnt vector will have sample.cut_num+2 elements
    int n_shells = sample.cut_num+2;
    
    //Initialize the shells_cnt vector
    shells_cnt.assign(n_shells, vector<int>());
    //hout << "n_shells="<<n_shells<<endl;
    
    //The coordinates of the middle point of the sample, boundary layer (outer shell), core shell, and
    //half step of observation window increment (or decrement) is needed several times
    //To reduce computations, these values are calcualted here and stored in arrays that are used
    //in the rest of the functions
    double midpoints[] = {
        sample.sample.poi_min.x+sample.sample.len_x/2.0,
        sample.sample.poi_min.y+sample.sample.wid_y/2.0,
        sample.sample.poi_min.z+sample.sample.hei_z/2.0};
    double boundary_layer[] = {
        sample.sample.poi_min.x+(sample.sample.len_x-sample.win_max_x)/2.0,
        sample.sample.poi_min.y+(sample.sample.wid_y-sample.win_max_y)/2.0,
        sample.sample.poi_min.z+(sample.sample.hei_z-sample.win_max_z)/2.0};
    double core[] = {
        sample.sample.poi_min.x+(sample.sample.len_x-sample.win_min_x)/2.0,
        sample.sample.poi_min.y+(sample.sample.wid_y-sample.win_min_y)/2.0,
        sample.sample.poi_min.z+(sample.sample.hei_z-sample.win_min_z)/2.0};
    double half_step[] = {sample.win_delt_x/2.0, sample.win_delt_y/2.0, sample.win_delt_z/2.0};
    
    
    //Scan all CNT points to determine which shells each CNTs occupies
    for (long int i = 0; i < (long int)points_cnt.size(); i++) {
        if (!Add_to_cnt_shells(midpoints, boundary_layer, core, half_step, points_cnt[i], n_shells, shells_cnt)) {
            hout << "Error in Generate_shells_and_structure when calling Add_to_cnt_shells."<< endl;
            return 0;
        }
    }
    
    //Scan all GNPs to determine which shells each CNTs occupies
    if (!Add_to_gnp_shells(midpoints, boundary_layer, core, half_step, gnps, n_shells, shells_gnps)) {
        hout << "Error in Generate_shells_and_structure when calling Add_to_gnp_shells."<< endl;
        return 0;
    }
    
    return 1;
}
//This function finds shell sub-region where a point is located.
//Then it adds the CNT number to the corresponding shell in the 2D vector shells_cnt
int Shells::Add_to_cnt_shells(const double midpoints[], const double boundary_layer[], const double core[], const double half_step[], const Point_3D &point, const int &n_shells, vector<vector<int> > &shells_cnt)
{
    
    //Find the shell that corresponds to the point
    int shell = Find_minimum_shell(midpoints, boundary_layer, core, half_step, point, n_shells);
            
    //Finally add the CNT on the corresponding sectioned domain
    //hout << "shell="<<shell<<endl;
    //Add the CNT number to the shell sub-region, only if:
    //it is empty
    //OR
    //the last CNT is not the current CNT (i.e., the CNT was already added to 'shell')
    if ( !shells_cnt[shell].size() || shells_cnt[shell].back() != point.flag ){
            shells_cnt[shell].push_back(point.flag);
    }
    
    return 1;
}
//This function finds the shell to which one point belongs to
//So it uses three times the function that finds the shell to which one coordinate belongs to
int Shells::Find_minimum_shell(const double midpoints[], const double boundary_layer[], const double core[], const double half_step[], const Point_3D &point, const int &n_shells)
{
    //Find the shell based on the x coordinate
    //hout << "x_in=";
    int shell_x = Find_shell(point.x, midpoints[0], boundary_layer[0], core[0], half_step[0], n_shells);
    //Find the shell based on the y coordinate
    //hout << "y_in=";
    int shell_y = Find_shell(point.y, midpoints[1], boundary_layer[1], core[1], half_step[1], n_shells);
    //Find the shell based on the z coordinate
    //hout << "z_in=";
    int shell_z = Find_shell(point.z, midpoints[2], boundary_layer[2], core[2], half_step[2], n_shells);
    
    //A CNT belongs to the outer-most shell from the three coordinates, that is,
    //the minimum shell number overall
    int shell = min(shell_z, min(shell_x, shell_y));
    
    return shell;
}
//This function finds the shell to which one coordinate belongs to
int Shells::Find_shell(const double &x_in, const double &x_m, const double &x_layer, const double &x_core, const double &dx_half, const int &n_shells)
{
    
    //Effective coordinate
    double x;
    //If the point is greater than the middle point, map it to the mirrored range
    if (x_in > x_m) x = 2*x_m - x_in;
    else x = x_in;
    
    //Check if x is in the outer shell, i.e., the boundary layer (it is in the boundary layer
    //when the point coordinate 'x' is below x_layer)
    //Also consider the points at the boundary layer to be part of the shell
    //hout<<x_in<<" x="<<x<<" shell";
    if (x < x_layer || abs(x - x_layer) < Zero) {
        //hout<<"m=0"<<endl;
        return 0;
    }
    //Check if x is in the inner shell
    else if( x > x_core ) {
        //hout<<"M="<<shells_cnt.size()-1<<endl;
        return n_shells-1;
    } else {
        //I need the integer part of (x - x_min)/dx_half + 1
        //The Zero is for floating point errors
        int shell = (int) ((x - x_layer + Zero)/dx_half) + 1;
        //hout<<"=f("<<((x - x_layer + Zero)/dx_half)+1<<")="<<shell<<endl;
        return shell;
    }
}
//This function finds the shell sub-regions that a GNP occupies.
int Shells::Add_to_gnp_shells(const double midpoints[], const double boundary_layer[], const double core[], const double half_step[], const vector<GNP> &gnps, const int &n_shells, vector<Shell> &shells_gnp)
{
    //Iterate over all GNPs
    for (int j = 0; j < (int)gnps.size(); j++) {
        
        //Variables to store the minimum and maximum shell numbers that the GNP occupies
        //Minimum and maximum shell numbers are initialized with the shell of the first GNP vertex
        shells_gnp[j].shell_min = Find_minimum_shell(midpoints, boundary_layer, core, half_step, gnps[j].vertices[0], n_shells);
        shells_gnp[j].shell_max = shells_gnp[j].shell_min;
        
        //Iterate over the remaining vertices of the GNP
        for (int i = 1; i < 8; i++) {
            
            //Find the shell that corresponds to vertex i
            int shell = Find_minimum_shell(midpoints, boundary_layer, core, half_step, gnps[j].vertices[i], n_shells);
            
            //Check if shell numbers need to be updated
            if (shell < shells_gnp[j].shell_min) {
                shells_gnp[j].shell_min = shell;
            }
            if (shell > shells_gnp[j].shell_max) {
                shells_gnp[j].shell_max = shell;
            }
        }
    }
    
    
    return 1;
}
