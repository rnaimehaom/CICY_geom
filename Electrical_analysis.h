//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Find the backbone and calculate the electrical resistivity and resistance on each direction
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef ELECTRICAL_ANALYSIS_h
#define ELECTRICAL_ANALYSIS_h

#include "Input_Reader.h"
#include "Backbone_Network.h"
#include "Cutoff_Wins.h"
#include "Direct_Electrifying.h"
#include "Hoshen_Kopelman.h"
#include "Printer.h"

//-------------------------------------------------------
class Electrical_analysis
{
public:
    
    //Constructor
    Electrical_analysis(){};
    
    int Perform_analysis_on_clusters(const int &iter, const cuboid &window, const Simu_para &simu_param, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const Visualization_flags &vis_flags, const Output_data_flags &out_flags, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps);
    int Get_number_of_clusters(const vector<vector<int> >& clusters_cnt, const vector<vector<int> >& clusters_gnp);
    int Clear_triangulations(const int& iter, const vector<vector<int> >& clusters_gnp, vector<GNP>& gnps);
    int Export_triangulations(const int &iter, const vector<int> &cluster_gnp, const vector<GNP> &gnps, const vector<Point_3D> &points_gnp);
    int Clear_triangulations_of_cluster(const vector<int>& cluster_gnp, vector<GNP>& gnps);
    int Electrical_resistance_along_each_percolated_direction (const int &R_flag, const int &n_cluster, const cuboid &window, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const Simu_para &simu_param, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps, vector<vector<double> > &paralel_resistors);
    int Vector_of_directions(const int &family, const Simu_para &simu_param, vector<int> &directions);
    int Calculate_parallel_resistor(const int &direction, const int &n_cluster, const Electric_para &electric_param, Direct_Electrifying *DEA, const vector<Point_3D> &points_cnt, const vector<double> &radii, vector<set<long int> > &elements, const vector<vector<int> > &clusters_cnt, const vector<vector<int> > &boundary_cnt, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const vector<vector<int> > &clusters_gnp, const vector<vector<int> > &boundary_gnp, vector<vector<double> > &paralel_resistors);
    int Get_boundaries_from_direction(const int &direction, int &b1, int &b2);
    int Currents_through_boundary_cnts(const Electric_para &electric_param, Direct_Electrifying *DEA, const vector<Point_3D> &points_cnt, const vector<double> &radii, vector<set<long int> > &elements, const vector<int> &boundary_cnt, double &I);
    int Current_of_element_in_boundary(const long int &P1, const long int &P2, const double &radius, Direct_Electrifying *DEA, const Electric_para &electric_param, const vector<Point_3D> &points_cnt, double &I);
    int Currents_through_boundary_gnps(const long int &node, const Electric_para &electric_param, Direct_Electrifying *DEA, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const vector<int> &boundary_gnp, double &I);
    int Calculate_percolated_families_fractions(const int &cnt_gnp_flag, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<GNP> &gnps, Hoshen_Kopelman *HoKo, Backbone_Network *BN);
    int Calculate_volume_of_non_percolated_cnts(const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<int> > &isolated_cnts, Backbone_Network *BN, double &np_cnts);
    int Calculate_volume_of_non_percolated_gnps(const vector<GNP> &gnps, const vector<vector<int> > &isolated_gnps, double &np_gnps);
    int Export_isolated_particles(const int &iter, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<vector<int> > &isolated_cnts, const vector<GNP> &gnps, const vector<vector<int> > &isolated_gnps);
    int Calculate_resistances_and_resistivities(const cuboid &window, const Electric_para &electric_param, const vector<vector<double> > &paralel_resistors);
    int Calculate_matrix_resistance(const int &direction, const cuboid &window, const double &matrix_resistivity, double &R_M);
    int Calculate_resistivity(const int &direction, const cuboid &window, const double &resistance, double &rho);
    

private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
