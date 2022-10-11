//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Read input data and save it into data structures
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef INPUTREADER_H
#define INPUTREADER_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Geometry_3D.h"

const double PI = 3.1415926535897932384626433832795;

//---------------------------------------------------------------------------
//Name of application case
struct App_name{
    string keywords;
    bool mark;
    string str;
};
//Name of simulation
struct Simu_para{
    string keywords;
    bool mark;
    string simu_name;
    //String to define whether a network is read from file, 
    //generated randomly or from seeds
    string create_read_network;
    //In case file is read from file (withput reading displacements from Abaqus)
    //this string specifies the type of file to be read: csv or dat (binary)
    string file_type;
    //Path to obd file or name if it is in the same folder as the executable
    string odb_file;
    //Name of the simulation step in Abaqus
    string step_name;
    //Number of frames in the Abaqus simulation (if reading displacements dirrectly
    //from a bynary file)
    int n_frames;
    //Flag to determine if a simplex of size three might be allowed to happen
    //when idintifying junction points on GNPs
    bool allow_simplex3;
    int sample_num;
    int simulation_scope;
    int resistances[3];
    int penetration_model_flag;
    string particle_type;
    //Seeds (positive integers) for the random number generators to generate a CNT and GNP networks
    vector<unsigned int> CNT_seeds, GNP_seeds;
    //Criterion to measure nanotube content: vol, wt
    string criterion;
    //Variables to store the volume and weight fractions
    double volume_fraction, weight_fraction;
    //Criterion for relationship of CNT to GNP content: mass_ratio, volume_ratio, 
    //relative volume fraction, relative weight fraction, density
    string mixed;
    //Mass ratio for mixed or hybrid particles
    double mass_ratio;
    //Mass ratio for mixed or hybrid particles
    double volume_ratio;
    //Relative volume fraction for mixed particles
    double relative_vol;
    //Relative weight fraction for mixed particles
    double relative_wt;
    //CNT density on GNPs
    double cnt_gnp_density;
    //Tolerance for error reduction in conjugate gradient
    double tolerance;
    //Maximum number of iterations for generating a CNT point
    int MAX_ATTEMPTS_CNT;
    //Maximum number of iterations for relocating a GNP
    int MAX_ATTEMPTS_GNP;
};
//The geometry of the sample
struct Geom_sample{
    string keywords;
    bool mark;
    double volume;
    double matrix_density;
    //Size for background grids (overlapping sub-regions)
    double gs_minx, gs_miny, gs_minz;
    //Overlapping of background grids
    double gs_overlap_cnt, gs_overlap_gnp;
    //Minimum, maximum, and step decrement of the observation window
    double win_min_x, win_min_y, win_min_z;
    double win_max_x, win_max_y, win_max_z;
    double win_delt_x, win_delt_y, win_delt_z;
    //Define the number of cutoff times (0: the maxmum size, n: the maxmum size - n*step_length(delta), n>=1)
    int cut_num;
    //Cuboid for the sample
    cuboid sample;
    //Cuboid for the extended domain for GNPs
    cuboid ex_dom_gnp;
    //Cuboid for the extended domain for CNTs
    cuboid ex_dom_cnt;
    //Cuboid for the non-penetration domain of CNTs
    cuboid np_domain;
    
    //Update the geometry of the observation window
    cuboid update_observation_window_geometry(const int &window, const string particle_type)
    {
        //Cuboid to be returned
        cuboid cub;
        
        //Update dimensions of the current observation window
        cub.len_x = win_max_x - ((double)window)*win_delt_x;
        cub.wid_y = win_max_y - ((double)window)*win_delt_y;
        //Set the dimension along z depending on the particle type
        cub.hei_z = (particle_type != "CNT_deposit")? win_max_z - ((double)window)*win_delt_z : sample.hei_z;
        //hout<<"len_x="<<cub.len_x<<" wid_y="<<cub.wid_y<<" hei_z="<<cub.hei_z<<endl;
        
        //These variables are the coordinates of the lower corner of the observation window
        cub.poi_min.x = sample.poi_min.x + (sample.len_x - cub.len_x)/2.0;
        cub.poi_min.y = sample.poi_min.y + (sample.wid_y - cub.wid_y)/2.0;
        cub.poi_min.z = sample.poi_min.z + (sample.hei_z - cub.hei_z)/2.0;
        //hout<<"poi_min="<<cub.poi_min.str()<<endl;
        
        //Boundaries with maximum values of coordinates
        cub.max_x = cub.poi_min.x + cub.len_x;
        cub.max_y = cub.poi_min.y + cub.wid_y;
        cub.max_z = cub.poi_min.z + cub.hei_z;
        //hout<<"max_x="<<cub.max_x<<" max_y="<<cub.max_y<<" max_z="<<cub.max_z<<endl;
        
        return cub;
    }
    
};
//The nanotube parameters in a network
struct Nanotube_Geo{
    string keywords;
    bool mark;
    //Criterion to measure nanotube content: vol, wt
    string criterion;
    //Initial growth direction type (random or specific)
    string dir_distrib_type;
    //Distribution type (uniform or normal) of the nanotube length
    string len_distrib_type;
    //Distribution type (uniform or normal) of the nanotube radius
    string rad_distrib_type;
    //Minimum nanotube length to keep a CNT close the the boundary (of the sample or observation window)
    string min_length_type;
    //Number of points to keep a CNT close the the boundary (of the sample or observation window)
    int min_points;
    //Length of CNT segment between two CNT consecutive points
    double step_length;
    //Angles for initial direction (used for 'specific' growth direction)
    double ini_theta, ini_phi;
    //Limits for the normal distribution of angle 'omega' in range [omega_a, omega_b]
    double omega_a, omega_b;
    //Length range (min, max) of nanotubes
    double len_min, len_max;
    //Radius range (min, max) of nanotubes
    double rad_min, rad_max;
    //Volume fraction of nanotubes
    double volume_fraction;
    //Actual volume of nanotubes
    double volume;
    //Weight fraction of nanotubes
    double weight_fraction;
    //Actual weight of nanotubes
    double weight;
    //Density of nanotubes
    double density;
};
//The nanotube parameters in a network
struct GNP_Geo{
    string keywords;
    bool mark;
    //Criterion to measure nanotube content: vol, wt
    string criterion;				
    //CNT growth type on GNP surfaces (parallel or independent)
    string growth_type;
    //Criterion for minimum GNP volume inside the sample
    string vol_in;
    //Minimum GNP volume inside the sample
    double min_vol_in;
    //GNP orientation type (random or specific)
    string orient_distrib_type;
    //Length distribution type (uniform or normal)
    string size_distrib_type;
    //Thickness distribution type (uniform or normal)
    string thick_distrib_type;
    //Angles for initial direction (used for 'specific' growth direction)
    double ini_theta, ini_phi;
    //GNP length range (min,max)
    double len_min, len_max;
    //GNP thickness range (min,max)
    double t_min, t_max;
    //CNT/GNP mass ratio
    double mass_ratio;
    //GNP volume fraction
    double volume_fraction;
    //Volume of GNPs inside the sample
    double volume;
    //GNP weight fraction
    double weight_fraction;
    //Weight of GNPs inside the sample
    double weight;
    //GNP density
    double density;
};
//Cutoff distances
struct Cutoff_dist{
    string keywords;
    bool mark;
    double van_der_Waals_dist;
    double tunneling_dist;
    int min_points;
    //Tolerance for considering two GNP junctions points (in the same GNP)
    //to be the same point
    double tol_gnp;
    //Squared value of tol_gnp
    //This is the quantity that will actually be used for comparison as squared
    //distances will be calculated to reduce computational time
    double tol_gnp2;
    //Maximum distance that a GNP is allowed to interpenetrate another GNP
    //(i.e., maximum distance to ignore interpenetrations)
    double PD_cutoff;
};
//Electrical parameters
struct Electric_para{
    string keywords;
    bool mark;
    //Applied voltage
    double applied_voltage;
    //CNT resistivity
    double resistivity_CNT;
    //GNP resistivities along the thicknes direction and along the surface
    double resistivity_GNP_t, resistivity_GNP_surf;
    //Resistivity of the polymer matrix
    double resistivity_matrix;
    //Type of junction resistance (constant or exponential)
    string junction_type;
    //Constants for tunneling
    double e_charge;
    //Constants using Li et al approach
    double e0_vacuum;
    double CNT_work_function;
    double K_polymer;
    //Constants using Hu et al approach
    double h_plank;
    double e_mass;
    double lambda_barrier;
    //Constants used to calculate the junction resistance when using the 
    //exponential equation by Hu et al.
    //This avoid to make operations that always give the same result
    double C1, C2;
    //Constant value for junction resistance
    double junction_resistance;
    //Scaling factor for reducing numerical errors
    //double scaling_R;
};
struct Visualization_flags{
    string keywords;
    bool mark;
    //Flag to export generated nanoparticles
    int generated_nanoparticles;
    //Flag to export clusters as obtained from the Hoshen-Kopelman algorithm
    int clusters;
    //Flag to export percolated clusters
    int percolated_clusters;
    //Flag to export the backbone
    int backbone;
    //Flag to export triangulations
    int triangulations;
    //Flags to export the sample and observation window domain
    int sample_domain;
    int window_domain;
};
struct Output_data_flags{
    string keywords;
    bool mark;
    //Flag to export separate files with fractions of CNTs and GNPs
    //(only when using mixed or hybrid particles)
    int cnt_gnp_flag;
    //Flag to export four vertices per GNP to generate a GNP network in Abaqus
    int gnp_data;
    //Precision (number of digits after the decimal point) used for
    //exporting the vertices
    int prec_gnp;
    //Flag to export CNT points
    int cnt_data;
    //Precision (number of digits after the decimal point) used for
    //exporting the CNT points
    int prec_cnt;
    
};
//---------------------------------------------------------------------------
class Input
{
public:
    //Data members
    App_name app_name;
    Simu_para simu_para;
    Geom_sample geom_sample;
    Nanotube_Geo nanotube_geo;
    Cutoff_dist cutoff_dist;
    Electric_para electric_para;
    GNP_Geo gnp_geo;
    Visualization_flags vis_flags;
    Output_data_flags out_flags;
    
    //Constructor
    Input(){};
    
    //Member functions
    int Data_Initialization();
    void Warning_message(const string &str);
    void Warning_message_already_input(const string &str);
    int Read_input_file(ifstream &infile);
    //Read the input data in a whole line (to skip over the comment line starting with a '%')
    string Get_Line(ifstream &infile)const;
private:
    //Member functions
    int Read_application(App_name &app_name, ifstream &infile);
    int Read_simulation_parameters(Simu_para &simu_para, ifstream &infile);
    int Read_sample_geometry(Geom_sample &geom_sample, ifstream &infile);
    int Read_cutoff_distances(Cutoff_dist &cutoff_dist, ifstream &infile);
    int Read_nanotube_geo_parameters(Nanotube_Geo &nanotube_geo, ifstream &infile);
    int Read_electrical_parameters(Electric_para &electric_para, ifstream &infile);
    int Read_gnp_geo_parameters(GNP_Geo &gnp_geo, ifstream &infile);
    int Read_visualization_flags(Visualization_flags &vis_flags, ifstream &infile);
    int Read_output_data_flags(Output_data_flags &out_flags, ifstream &infile);
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
