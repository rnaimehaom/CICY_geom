//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
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
    
    //vector<double> resistances;
    vector<vector<long int> > elements; //This vector will store the elements. It is needed to trim the CNTs
    vector<vector<long int> > elements_tunnel, elements_mixed_tunnel, elements_gnp_tunnel; //This vector will store the tunnel elements. It is needed to calculate the zero-current cutoff
    vector<int> LM_matrix;//Local mapping matrix. It maps from point number to node number. It is also used to calculate the currents
    vector<int> LM_matrix_gnp;//Local mapping matrix. It maps from point number to node number. It is also used to calculate the currents
    vector<vector<vector<int> > > boundary_node_map;
    int reserved_nodes; //This is the number of nodes assigned to the boudaries with prescribed voltage
    long int last_cnt_node; //This is the las CNT node number
    
    //Constructor
    Direct_Electrifying(){};
    
    //Member Functions
    int Compute_voltage_field(const int &n_cluster, const int &R_flag, const Electric_para &electric_param, const Cutoff_dist &cutoffs, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps);
    int Get_global_nodes(const int &family);
    int LM_matrix_for_cnts(const int &n_cluster, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, long int &global_nodes);
    int Map_points_at_boundaries(const int &family, const vector<vector<long int> > &boundary_pts, map<long int, long int> &LMM);
    int Get_vector_of_boundaries(const int &family, vector<int> &boundaries);
    int LM_matrix_for_gnps(const int &n_cluster, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure_gnp, long int &global_nodes);
    int Fill_sparse_stiffness_matrix(const int &R_flag, const long int &nodes, const long int &reserved_nodes, const double &d_vdw, const int &n_cluster, const Electric_para &electric_param, Hoshen_Kopelman *HoKo, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, vector<double> &P, vector<double> &R, vector<double> &VEF);
    int Fill_2d_matrices_cnts(const int &R_flag, const int &n_cluster, const Electric_para &electric_param,  Hoshen_Kopelman *HoKo, const vector<Point_3D> &points_cnt, const vector<double> &radii, vector<map<long int, double>> &col_values, vector<double> &diagonal);
    int Calculate_resistance_cnt(const int &R_flag, const vector<Point_3D> &points_cnt, const long int &P1, const long int &P2, const double &radius, const double &resistivity, double &Re);
    int Add_new_elements_to_2d_sparse_matrix(const long int &node1, const long int &node2, const double &Re, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int Add_to_existing_elements_in_2d_sparse_matrix(const long int &node1, const long int &node2, const double &Re, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int Fill_2d_matrices_cnt_junctions(const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_cnt_junctions_i, const vector<Junction> &junctions_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const map<long int, long int> &LMM_cnts, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int Calculate_junction_resistance(const Junction &j, const double &d_vdw, const double &rad1, const Point_3D &P1, const double &rad2, const Point_3D &P2, const struct Electric_para &electric_param, double &Re);
    int Fill_2d_matrices_mixed_junctions(const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_mix_junctions_i, const vector<Junction> &junctions_mixed, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const map<long int, long int> &LMM_cnts, const map<long int, long int> &LMM_gnps, vector<map<long int, double> > &col_values, vector<double> &diagonal, map<long int, double> &points_cnt_rad);
    int Fill_2d_matrices_gnp(const int &R_flag, const Electric_para &electric_param, const vector<int> &cluster_gnp, const vector<Point_3D> &points_gnp, const vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, const map<long int, long int> &LMM_gnps, map<long int, double> &points_cnt_rad, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int Calculate_resistance_gnp(const Point_3D &P1, const Point_3D &P2, const double &rad1, const double &rad2, const struct Electric_para &electric_param, double &Re);
    int Fill_2d_matrices_gnp_junctions(const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_gnp_junctions_i, const vector<Junction> &junctions_gnp, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const map<long int, long int> &LMM_gnps, vector<map<long int, double> > &col_values, vector<double> &diagonal);
    int From_2d_to_1d_vectors(const long int &reserved_nodes, const vector<map<long int, double> > &col_values, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal);
    int Set_up_residual_and_search_direction(const int &R_flag, const long int nodes, const long int &reserved_nodes, const Electric_para &electric_param, const vector<vector<double> > &KEFT, vector<double> &P, vector<double> &R, vector<double> &VEF);
    int Get_voltage_vector(const double &volts, const long int &reserved_nodes, vector<double> &VEF);
    int Solve_DEA_equations_CG_SSS(const long int &nodes, const long int &reserved_nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, vector<double> &P, vector<double> &R, vector<double> &VEF);
    void Jacobi_preconditioner(const vector<double> &diagonal, vector<double> &M_inv);
    void Apply_preconditioner(const vector<double> &M_inv, const vector<double> &R, vector<double> &P, vector<double> &Y);
    void spM_V_SSS(const vector<double> &V, const vector<long int> &rowptr, const vector<long int> &colind, const vector<double> &diagonal, const vector<double> &values, vector<double> &R);
    int V_dot_v(const vector<double> &A, const vector<double> &B, double &dot_product);
    void V_plus_aW(const vector<double> &W, const double &a, vector<double> &V);
    void W_plus_aV(const vector<double> &W, const double &a, vector<double> &V);
    void Componentwise_multiply(const vector<double> &vector_in1, const vector<double> &vector_in2, vector<double> &vector_out);
        
        
        
        
        
        
    
    int Calculate_voltage_field(const int &family, const int &n_cluster, const int &R_flag, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &point_list_gnp, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, vector<GCH> &hybrid_particles);
    void Initialize_boundary_node_map();
    int Get_LM_matrices(const int &family, const int &n_cluster, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure, const vector<vector<long int> > &structure_gnp, int &global_nodes, vector<int> &LM_matrix, vector<int> &LM_matrix_gnp, vector<vector<long int> > &elements, vector<vector<long int> > &gnp_triangulation_points);
    int Get_LM_matrix_cnts_only(const int &family, const int &n_cluster, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure, int &global_nodes, vector<int> &LM_matrix, vector<vector<long int> > &elements);
    int Fill_mixed_contact_flags(const int &family, const vector<contact_pair> &mixed_contacts, const vector<vector<short int> > &boundary_flags_cnt, vector<long int> &cnt_contacts_point, vector<long int> &gnp_contacts_point);
    int Fill_gnp_contact_flags(const vector<contact_pair> &gnp_contacts, vector<long int> &gnp_contacts_point);
    void Add_point_to_LM_matrix(long int P, int family, const vector<vector<short int> > &boundary_flags, int &global_nodes, vector<int> &LM_matrix);
    int Is_in_relevant_boundary(int family, short int boundary_node);
    int Get_boundary_node(const vector<short int> &boundary_flag, const int &family);
    int Fill_sparse_stiffness_matrix(const int &R_flag, const int &nodes, const double &d_vdw, const int &n_cluster, Hoshen_Kopelman *HoKo, const vector<vector<long int> > &structure, const vector<vector<long int> > &elements, const vector<int> &LM_matrix, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<vector<long int> > &gnp_triangulation_points, const vector<int> &LM_matrix_gnp, const vector<Point_3D> &point_list_gnp, vector<GCH> &hybrid_particles, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, const struct Electric_para &electric_param);
    int Fill_2d_matrices_cnts(const int &R_flag, const vector<vector<long int> > &elements, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &cluster, const vector<int> &LM_matrix, const struct Electric_para &electric_param, const double &d_vdw, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point);
    double Calculate_resistance_cnt(const vector<Point_3D> &point_list, const long int &P1, const long int &P2, const double &radius, const double &resistivity);
    void Add_elements_to_sparse_stiffness(const long int &node1, const long int &node2, const double &Re, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal);
    void Check_for_other_elements(const int &R_flag, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &LM_matrix, const struct Electric_para &electric_param, const double d_vdw, const long int &P1, const long int &node1, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point);
    double Calculate_resistance_tunnel(const double &rad1, const Point_3D &P1, const double &rad2, const Point_3D &P2, const struct Electric_para &electric_param, const double &d_vdw);
    double Calculate_resistance_tunnel(const int &tunnel_flag, const GCH &hybrid1, const Point_3D &P1, const double &rad2, const Point_3D &P2, const struct Electric_para &electric_param, const double &d_vdw);
    double Calculate_distance_tunnel_point_gnp(const GCH &hybrid1, const Point_3D &P1, const Point_3D &P2);
    double Exponential_tunnel(const struct Electric_para &electric_param, const double &d_vdw, const double &A, double &separation);
    void Remove_from_vector(long int num, vector<long int> &vec);
    int Fill_2d_matrices_gnp(const int &R_flag, const vector<int> &cluster_gnp, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &LM_matrix, const vector<Point_3D> &point_list_gnp, const vector<vector<long int> > &gnp_triangulation_points, const vector<int> &LM_matrix_gnp, const struct Electric_para &electric_param, vector<GCH> &hybrid_particles, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal);
    int Same_boundary_node_preprocessing(const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<int> &LM_matrix, const vector<int> &LM_matrix_gnp, const vector<long int> &gnp_triangulation_points, GCH &hybrid, vector<int> &tmp_cnts);
    int Only_add_values(const vector<int> &LM_matrix, const vector<short int> &edge_flags, const long int &node1, const int &particle1, const long int &node2, const int &particle2);
    double Calculate_resistance_gnp(const int &flag, const Point_3D &P1, const Point_3D &P2, const double &rad1, const double &rad2, const GCH &hybrid, const struct Electric_para &electric_param);
    int Flags_gnps_inside(const vector<int> &clusters_gnp, vector<short int> &gnps_inside_flags);
    int Fill_2d_matrices_gnp_gnp_tunnels(const int &R_flag, const vector<contact_pair> &gnp_contacts, const vector<Point_3D> &point_list_gnp, const vector<int> &LM_matrix_gnp, const vector<int> &clusters_gnp, vector<short int> &gnps_inside_flags, const struct Electric_para &electric_param, const double &d_vdw, vector<GCH> &hybrid_particles, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal);
    int Fill_2d_matrices_cnt_gnp_tunnels(const int &R_flag, const vector<contact_pair> &mixed_contacts, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &LM_matrix, const vector<Point_3D> &point_list_gnp, const vector<int> &LM_matrix_gnp, vector<short int> &gnps_inside_flags, const struct Electric_para &electric_param, const double &d_vdw, vector<GCH> &hybrid_particles, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal);
    int Check_repeated_col_ind_2d(const int &nodes, const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<Point_3D> &point_list, const vector<int> &LM_matrix, const vector<Point_3D> &point_list_gnp, const vector<int> &LM_matrix_gnp, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d);
    void Find_repeated_elements(const vector<long int> &vector_in, vector<long int> &elements, vector<int> &indices);
    int Export_matlab_sparse_matrix(const vector<vector<long int> > &col_ind_2d, const vector<vector<double> > &values_2d, const vector<double> &diagonal, const string &filename);
    void From_2d_to_1d_vectors(vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal);int Solve_DEA_equations_CG_SSS(const int &R_flag, long int nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, const struct Electric_para &electric_param, vector<vector<double> > &KEFT);
    void Get_voltage_vector(const double &nodes, vector<double> &voltages);
    void Conjugate_gradient(long int nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, vector<double> &R, vector<double> &P);
    double V_dot_v(const vector<double> &A, const vector<double> &B);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
