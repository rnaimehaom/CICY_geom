//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
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
    //Data Member
    vector<double> resistors;
    
    int Perform_analysis_on_clusters(const int &avoid_resistance_flag, const int &vtk_flag, const cuboid &window, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const Visualization_flags &vis_flags, const vector<int> &family, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices, vector<vector<int> > &gnp_dead_indices, vector<vector<int> > &gnp_indices);
    int Electrical_resistance_along_each_percolated_direction (const int &R_flag, const int &n_cluster, const int &family, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps, vector<vector<double> > &paralel_resistors);
    int Vector_of_directions(const int &family, vector<int> &directions);
    int Calculate_parallel_resistor(const int &direction, const int &n_cluster, const Electric_para &electric_param, Direct_Electrifying *DEA, const vector<Point_3D> &points_cnt, const vector<double> &radii, vector<set<long int> > &elements, const vector<vector<int> > &clusters_cnt, const vector<vector<int> > &boundary_cnt, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const vector<vector<int> > &clusters_gnp, const vector<vector<int> > &boundary_gnp, vector<vector<double> > &paralel_resistors);
    int Get_boundaries_from_direction(const int &direction, int &b1, int &b2);
    int Currents_through_boundary_cnts(const Electric_para &electric_param, Direct_Electrifying *DEA, const vector<Point_3D> &points_cnt, const vector<double> &radii, vector<set<long int> > &elements, const vector<int> &boundary_cnt, double &I);
    int Current_of_element_in_boundary(const long int &P1, const long int &P2, const double &radius, Direct_Electrifying *DEA, const Electric_para &electric_param, const vector<Point_3D> &points_cnt, double &I);
    int Currents_through_boundary_gnps(const int &node, const Electric_para &electric_param, Direct_Electrifying *DEA, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const vector<int> &boundary_gnp, double &I);
    int Calculate_percolated_families_fractions(const int &cnt_gnp_flag, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<GNP> &gnps, Hoshen_Kopelman *HoKo, Backbone_Network *BN);
    int Calculate_volume_of_non_percolated_cnts(const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<int> > &isolated_cnts, Backbone_Network *BN, double &np_cnts);
    int Calculate_volume_of_non_percolated_gnps(const vector<GNP> &gnps, const vector<vector<int> > &isolated_gnps, double &np_gnps);
    int Calculate_resistances_and_resistivities(const cuboid &window, const Electric_para &electric_param, const vector<vector<double> > &paralel_resistors);
    int Calculate_matrix_resistance(const int &direction, const cuboid &window, const double &matrix_resistivity, double &R_M);
    int Calculate_resistivity(const int &direction, const cuboid &window, const double &resistance, double &rho);
    
    
    
    
    
    
    
    //Deprecated:
    int Perform_analysis_on_clusters(const int &avoid_resistance_flag, const vector<int> &family, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &point_list_gnp, const struct Geom_sample &window, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, vector<GCH> &hybrid_particles, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices, vector<vector<int> > &gnp_dead_indices, vector<vector<int> > &gnp_indices);
    int Electrical_analysis_along_each_percolated_direction (const int &R_flag, const int &n_cluser, const vector<int> &family, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, Backbone_Network *Backbonet, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<long int> > &contacts_initial, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &point_list_gnp, vector<GCH> &hybrid_particles, vector<vector<double> > &paralel_resistors);
    int Convert_index_to_structure(const vector<int> &cluster, const vector<vector<long int> > &indices, vector<vector<long int> > &structure, vector<int> &backbone_cnts);
    int Update_hybrids(const vector<int> &cluster_gch, const vector<vector<long int> > &structure, const vector<vector<long int> > &backbone_structure, vector<GCH> &hybrid_particles);
    int Update_vectors_for_hoko_cutwins(const int &n_cnts, const int &n_gnps, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, Hoshen_Kopelman *HoKo_Re, Cutoff_Wins *Cutwins_Re, const vector<vector<long int> > &contacts_initial, const vector<vector<long int> > &structure, const vector<int> &backbone_cnts, const vector<int> &percolated_gnps);
    int Calculate_parallel_resistor(const int &direction, Direct_Electrifying * DEA, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<int> > &clusters_cnt, const vector<vector<int> > &boundary_cnt, const vector<Point_3D> &point_list_gnp, const vector<GCH> &hybrid_particles, const vector<vector<int> > &clusters_gnp, const vector<vector<int> > &boundary_gnp, const vector<vector<short int> > &boundary_flags_gnp, const struct Electric_para &electric_param, vector<vector<double> > &paralel_resistors);
    double Current_of_element_in_boundary(const long int &P1, const long int &P2, const double &radius, Direct_Electrifying *DEA, const struct Electric_para &electric_param, const vector<Point_3D> &point_list);
    double Current_of_edges_in_boundary(const int &side, Direct_Electrifying *DEA, const struct Electric_para &electric_param, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<Point_3D> &point_list_gnp, const GCH &hybrid);
    int Calculate_matrix_resistances(const double &matrix_resistivity, const struct Geom_sample &window, vector<double> &matrix_resistances);
    int Calculate_resistances_and_resistivities(const struct Geom_sample &window, const vector<double> &matrix_resistances, const vector<vector<double> > &paralel_resistors, vector<double> &resistors, vector<double> &resistivities);

private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
