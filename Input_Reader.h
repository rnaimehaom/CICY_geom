//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
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
    int penetration_model_flag;
    string particle_type;
    //Seeds (positive integers) for the random number generators to generate a CNT and GNP networks
    vector<unsigned int> CNT_seeds, GNP_seeds;
    //Criterion to measure nanotube content: vol, wt
    string criterion;
    //Variables to store the volume and weight fractions
    double volume_fraction, weight_fraction;
    //Criterion for relationship of CNT to GNP content: mr, dens
    string mixed;
    //Mass ratio for mixed or hybrid particles
    double mass_ratio;
    //Mass ratio for mixed or hybrid particles
    double volume_ratio;
    //CNT density on GNPs
    double cnt_gnp_densinty;
};
//The geometry of the sample
struct Geom_sample{
    string keywords;
    bool mark;
    //Lower left corner of sample
    Point_3D origin;
    //Length, width and height of sample domain
    double len_x, wid_y, hei_z;
    //Coordinates of the sample's boundaries opposite to those given by the coordinates of origin
    double x_max, y_max, z_max;
    double volume;
    double matrix_density;
    //Minimum size for background grids (overlapping sub-regions)
    double gs_minx, gs_miny, gs_minz;
    //Overlapping of background gids
    double gs_overlap;
    //Minimum, maximum, and step decrement of the observation window
    double win_min_x, win_min_y, win_min_z;
    double win_max_x, win_max_y, win_max_z;
    double win_delt_x, win_delt_y, win_delt_z;
    //Define the number of cutoff times (0: the maxmum size, n: the maxmum size - n*step_length(delta), n>=1)
    int cut_num;
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
    //Step length (unit: micromether) for discretization of the GNP
    double discr_step_length;
    //Angles for initial direction (used for 'specific' growth direction)
    double ini_theta, ini_phi;
    //GNP length range (min,max)
    double len_min, len_max;
    //GNP thickness range (min,max)
    double t_min, t_max;
    //CNT/GNP mass ratio
    double mass_ratio;
    //GNp volume fraction
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
//The parameters of nanotube agglomerates
struct Agglomerate_Geo{
    string keywords;
    int print_key;								//0 denotes "no print"; 1 denotes "only print the nanotubes in the ellipsoids"; 2 denotes "print the nanotubes in the ellipsolds and the surfaces of all ellipsoids".
    bool mark;
    double volf_clust;						//Define the volume fraction of nanotubes in clusters
    double vol_fra_criterion;			//Define the volume fraction of clusters in the sample
    double amin;								//Define the minimum value of long axis of a cluster ellipsoid
    double amax;							//Define the maximum value of long axis of a cluster ellipsoid
    double bmin;								//Define the minimum value of middle axis of a cluster ellipsoid
    double cmin;								//Define the minimum value of short axis of a cluster ellipsoid
    double growth_probability;		//Define the growth probability of nanotubes in a cluster
    double cnt_real_volume;			//Define the real volume of nanotubes in clusters
    vector<struct elliparam> ellips;  //Define the vector of ellipsoids for nanotube cluster zones
};
//The cutoff distances
struct Cutoff_dist{
    string keywords;
    bool mark;
    double van_der_Waals_dist;
    double tunneling_dist;
};
//The electrical parameters
struct Electric_para{
    string keywords;
    bool mark;
    double applied_voltage;			//Define the magnitude of the applied voltage
    double resistivity_CNT;			//Define the resistivity value of the carbon fiber
    double resistivity_GNP_t, resistivity_GNP_surf; //Define the resistivity value of the GNP along the thicknes direction and along the surface
    double resistivity_matrix;      //Define the resistivity of the polymer matrix
    string junction_type;           //Define the type of junctino resistance
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
    double junction_resistance;     //Define the constant value for junction resistance
};
struct Tecplot_flags{
    string keywords;
    bool mark;
    int generated_cnts;         //Flag to export generated CNTs
    int generated_gnps;         //Flag to export generated GNPs
    int clusters;               //Flag to export clusters as obtained from the Hoshen-Kopelman algorithm
    int percolated_clusters;    //Flag to export percolated clusters
    int backbone;               //Flag to export the backbone
    int triangulations;         //Flag to export triangulations
};
//---------------------------------------------------------------------------
class Input
{
public:
    //Data members
    struct App_name app_name;
    struct Simu_para simu_para;
    struct Geom_sample geom_sample;
    struct Nanotube_Geo nanotube_geo;
    struct Agglomerate_Geo agg_geo;
    struct Cutoff_dist cutoff_dist;
    struct Electric_para electric_para;
    struct GNP_Geo gnp_geo;
    struct Tecplot_flags tec360_flags;
    
    //Constructor
    Input(){};
    
    //Member functions
    int Data_Initialization();								//Initialize data
    void Warning_message(const string &str);
    void Warning_message_already_input(const string &str);
    int Read_Infile(ifstream &infile);				//Read data
    string Get_Line(ifstream &infile)const;		//Read the input data in a whole line (to skip over the comment line starting with a '%')
private:
    //Member functions
    int Read_application(struct App_name &app_name, ifstream &infile);
    int Read_simulation_parameters(struct Simu_para &simu_para, ifstream &infile);
    int Read_sample_geometry(struct Geom_sample &geom_sample, ifstream &infile);
    int Read_cutoff_distances(struct Cutoff_dist &cutoff_dist, ifstream &infile);
    int Read_nanotube_geo_parameters(struct Nanotube_Geo &nanotube_geo, ifstream &infile);
    int Read_agg_geo_parameters(struct Agglomerate_Geo &agg_geo, ifstream &infile);
    int Read_electrical_parameters(struct Electric_para &electric_para, ifstream &infile);
    int Read_gnp_geo_parameters(struct GNP_Geo &gnp_geo, ifstream &infile);
    int Read_tecplot_flags(struct Tecplot_flags &tec360_flags, ifstream &infile);
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
