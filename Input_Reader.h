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
const int MAX_ATTEMPTS = 5;

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
    string create_read_network;
    int sample_num;
    int avoid_resistance;
    int resistances[3];
    int penetration_model_flag;
    string particle_type;
    //Seeds (positive integers) for the random number generators to generate a CNT and GNP networks
    vector<unsigned int> CNT_seeds, GNP_seeds;
    //Criterion to measure nanotube content: vol, wt
    string criterion;
    //Variables to store the volume and weight fractions
    double volume_fraction, weight_fraction;
    //Criterion for relationship of CNT to GNP content: mass_ratio, volume_ratio, density
    string mixed;
    //Mass ratio for mixed or hybrid particles
    double mass_ratio;
    //Mass ratio for mixed or hybrid particles
    double volume_ratio;
    //CNT density on GNPs
    double cnt_gnp_density;
    //Tolerance for error reduction in conjugate gradient
    double tolerance;
};
//The geometry of the sample
struct Geom_sample{
    string keywords;
    bool mark;
    double volume;
    double matrix_density;
    //Minimum size for background grids (overlapping sub-regions)
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
    //Angle 'omega' for the normal distribution range [-omega, omega] of growth direction
    double angle_max;
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
    //Constant value for junction resistance
    double junction_resistance;
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
    int gnp_4p;
    //Precision (number of digits after the decimal point) used for
    //exporting the vertices
    int prec_gnp;
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
