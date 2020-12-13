//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Determine the backbone network and dead branches in the percolation network
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef BACKBONE_NETWORK_H
#define BACKBONE_NETWORK_H

#include "Direct_Electrifying.h"
#include "Hoshen_Kopelman.h"
#include "Input_Reader.h"
#include "Printer.h"
#include "VTK_Export.h"
#include <algorithm>
#include <set>
#include <map>
#include <iterator>


//-------------------------------------------------------
class Backbone_Network
{
public:
    //Variables to store the volumes of particles
    vector<double> volumes_cnt;
    vector<double> volumes_gnp;
    //Variable to store the volumes of dead particles
    vector<double> dead_branches;
    vector<double> dead_gnps;
    //Variables to store total volumes per family
    vector<double> total_volumes;
    vector<double> dead_particles;
    
    //Constructor
    Backbone_Network(){
        
        //Initialize vectors
        volumes_cnt.assign(8,0);
        volumes_gnp.assign(8,0);
        dead_branches.assign(7,0);
        dead_gnps.assign(7,0);
        total_volumes.assign(8,0);
        dead_particles.assign(7,0);
    };
    
    //Member Functions
    int Determine_backbone_network(const int &n_cluster, const int &R_flag, const int &avoid_resistance_flag, const int &vtk_flag, const vector<double> &voltages, const map<long int, long int> &LMM_cnts, const map<long int, long int> &LMM_gnps, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, Hoshen_Kopelman *HoKo);
    int Find_zero_current(const int &n_cluster, const int &R_flag, const vector<double> &voltages, const map<long int, long int> &LMM_cnts, const map<long int, long int> &LMM_gnps, Hoshen_Kopelman *HoKo, const Cutoff_dist &cutoffs, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, vector<vector<double> > &currents_cnt, vector<vector<double> > &currents_gnp, double &zero_current);
    int Zero_current_in_same_particle_junctions_unit_resistor(const vector<int> &junction_cluster, const vector<Junction> &junctions, const vector<double> &voltages, const map<long int, long int> &LMM, double &I_max);
    int Zero_current_in_same_particle_junctions_unit_resistor_test(const vector<int> &junction_cluster, const vector<Junction> &junctions, const vector<double> &voltages, const map<long int, long int> &LMM, double &I_max, vector<double> &currents);
    int Find_backbone_and_fractions_cnts(const int &n_cluster, const int &avoid_resistance_flag, const int &vtk_flag, const double &zero_current, const vector<vector<double> > &currents_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, Hoshen_Kopelman *HoKo);
    int Calculate_cnt_volumes(const int &n_cluster, const int &avoid_resistance_flag, const int &vtk_flag, const int &CNTi, const int &idx1, const int &idx2, const vector<Point_3D> &points_cnt, const double &radius, Hoshen_Kopelman *HoKo, vector<long int> &dead_branches_i, vector<long int> &backbone_i);
    int CNT_volume_between_two_points(const long int &P1, const long int P2, const double &radius, const vector<Point_3D> &points_cnt, double &volume);
    int Find_backbone_and_fractions_gnps(const int &n_cluster, const int &avoid_resistance_flag, const int &vtk_flag, const double &zero_current, const vector<vector<double> > &currents_gnp, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, Hoshen_Kopelman *HoKo);
    int Remove_junctions_with_dead_particles_cnts(const int &n_cluster, Hoshen_Kopelman *HoKo);
    int Remove_junctions_with_dead_particles_gnps(const int &n_cluster, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, Hoshen_Kopelman *HoKo);
    int Remove_point_from_vector(const long int &P, vector<long int> &structure_i);
    int Remove_junctions_with_dead_particles_mixed(const int &n_cluster, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, Hoshen_Kopelman *HoKo);
    
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
