//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Implementation of the Direct Electrifying Algorithm
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef DIRECTELECTRIFYING_H
#define DIRECTELECTRIFYING_H

#include "Cutoff_Wins.h"
#include "Hoshen_Kopelman.h"
#include "Input_Reader.h"
#include "Printer.h"
#include "Triangulation.h"
#include <algorithm>

//-------------------------------------------------------
class Direct_Electrifying
{
public:
    //Data Member
    vector<double> voltages;
    //Variables for the local mapping matrices (LMM)
    //These matricew map from point number to node number in the stiffness matrix
    map<long int, long int> LMM_cnts;
    map<long int, long int> LMM_gnps;
    //Set used to determine the nodes from CNT points in mixed junctions
    map<long int, double> points_cnt_rad;

    //Vector to output all resistors in the stiffness matrix
    //vector<double> all_resistors;
    
    //Constructor
    Direct_Electrifying(){};
    
    //Member Functions
    int Compute_voltage_field(const int &n_cluster, const int &R_flag, const cuboid &window, const Simu_para &simu_para, const Electric_para &electric_param, const Cutoff_dist &cutoffs, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps);
    int Get_global_nodes(const int &family);
    int LM_matrix_for_cnts(const int &n_cluster, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, long int &global_nodes);
    int Map_points_at_boundaries(const int &family, const vector<vector<long int> > &boundary_pts, map<long int, long int> &LMM);
    int Get_vector_of_boundaries(const int &family, vector<int> &boundaries);
    int LM_matrix_for_gnps(const int &n_cluster, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure_gnp, long int &global_nodes);
    int Fill_sparse_stiffness_matrix(const int &R_flag, const long int &nodes, const long int &reserved_nodes, const double &d_vdw, const int &n_cluster, const Electric_para &electric_param, Hoshen_Kopelman *HoKo, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, vector<double> &R, vector<double> &VEF);
    int Fill_2d_matrices_cnts(const int &R_flag, const int &n_cluster, const Electric_para &electric_param,  Hoshen_Kopelman *HoKo, const vector<Point_3D> &points_cnt, const vector<double> &radii, vector<map<long int, double>> &col_values, vector<double> &diagonal);
    int Three_points_in_element(const int &R_flag, const Electric_para &electric_param, const set<long int> &cnt_element, const vector<Point_3D> &points_cnt, const double &radius, vector<map<long int, double>> &col_values, vector<double> &diagonal);
    int Calculate_resistance_cnt(const int &R_flag, const vector<Point_3D> &points_cnt, const long int &P1, const long int &P2, const double &radius, const double &resistivity, double &Re);
    int Add_new_elements_to_2d_sparse_matrix(const long int &node1, const long int &node2, const double &Re_inv, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int Add_to_existing_elements_in_2d_sparse_matrix(const long int &node1, const long int &node2, const double &Re_inv, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int Add_new_or_to_existing_elements_in_2d_sparse_matrix(const long int &node1, const long int &node2, const double &Re_inv, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int Fill_2d_matrices_cnt_junctions(const int &R_flag, const long int &reserved_nodes, const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_cnt_junctions_i, const vector<Junction> &junctions_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const map<long int, long int> &LMM_cnts, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int Calculate_junction_resistance(const Junction &j, const double &d_vdw, const double &l1, const Point_3D &P1, const double &l2, const Point_3D &P2, const struct Electric_para &electric_param, double &Re);
    int Fill_2d_matrices_mixed_junctions(const int &R_flag, const long int &reserved_nodes, const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_mix_junctions_i, const vector<Junction> &junctions_mixed, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const map<long int, long int> &LMM_cnts, const map<long int, long int> &LMM_gnps, vector<map<long int, double> > &col_values, vector<double> &diagonal, map<long int, double> &points_cnt_rad);
    int Fill_2d_matrices_gnp(const int &R_flag, const Electric_para &electric_param, const vector<int> &cluster_gnp, const vector<Point_3D> &points_gnp, const vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, const map<long int, long int> &LMM_gnps, map<long int, double> &points_cnt_rad, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int Calculate_resistance_gnp(const Point_3D &P1, const Point_3D &P2, const double &rad1, const double &rad2, const struct Electric_para &electric_param, double &Re);
    int Fill_2d_matrices_gnp_junctions(const int &R_flag, const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_gnp_junctions_i, const vector<Junction> &junctions_gnp, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const map<long int, long int> &LMM_gnps, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int From_2d_to_1d_vectors(const long int &reserved_nodes, const long int &nodes, const int &R_flag, const Electric_para &electric_param, const vector<map<long int, double> > &col_values, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, vector<double> &R, vector<double> &VEF);
    int Initial_guess_for_CG(const int &n_cluster, const long int &reserved_nodes, const int &family, const cuboid &window_geom, const double &V_app, const vector<vector<int> > &clusters_cnt, const vector<set<long int> > &elements_cnt, const vector<Point_3D> &points_cnt, const vector<vector<int> > &clusters_gnp, const vector<Point_3D> &points_gnp, vector<double> &R, vector<double> &V_guess);
    double Select_coordinate(const double& family, const Point_3D& P);
    int Initial_guess_for_cnt_nodes(const int& n_cluster, const long int& reserved_nodes, const int& family, const double& V_P_coord, const double& max_coord, const vector<vector<int> >& clusters_cnt, const vector<set<long int> >& elements_cnt, const vector<Point_3D>& points_cnt, vector<double>& V_guess);
    int Get_voltage_vector(const long int &reserved_nodes, const long int &nodes, const int &R_flag, const Electric_para &electric_param, vector<double> &VEF);
    int Update_residual_vector(const vector<long int>& col_ind, const vector<long int>& row_ptr, const vector<double>& values, const vector<double>& diagonal, const vector<double>& V_guess, vector<double>& R);
    int Solve_DEA_equations_CG_SSS(const long int &nodes, const long int &reserved_nodes, const double &tolerance, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, vector<double> &R, vector<double> &VEF, vector<double> &voltages_sol);
    void Jacobi_preconditioner(const vector<double> &diagonal, vector<double> &M_inv);
    void Apply_preconditioner(const vector<double> &M_inv, const vector<double> &R, vector<double> &P, vector<double> &Y);
    void spM_V_SSS(const vector<double> &V, const vector<long int> &rowptr, const vector<long int> &colind, const vector<double> &diagonal, const vector<double> &values, vector<double> &R);
    int V_dot_v(const vector<double> &A, const vector<double> &B, double &dot_product);
    void V_plus_aW(const vector<double> &W, const double &a, vector<double> &V);
    void W_plus_aV(const vector<double> &W, const double &a, vector<double> &V);
    void Componentwise_multiply(const vector<double> &vector_in1, const vector<double> &vector_in2, vector<double> &vector_out);
        
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
