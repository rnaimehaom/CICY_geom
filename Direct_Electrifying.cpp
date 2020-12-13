//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Implementation of the Direct Electrifying Algorithm
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Direct_Electrifying.h"

//Calculate the voltage values at contact points and endpoints
int Direct_Electrifying::Compute_voltage_field(const int &n_cluster, const int &R_flag, const Electric_para &electric_param, const Cutoff_dist &cutoffs, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps)
{
    //First we need to prepare the matrices that the direct electrifying needs.
    //The number of reserved nodes is calculated.
    //These is the number of boundaries with prescribed voltage.
    //hout<<"Get_global_nodes"<<endl;
    int reserved_nodes = Get_global_nodes(HoKo->family[n_cluster]);
    if (reserved_nodes == -1) {
        hout<<"Error in Compute_voltage_field, invalid family: "<<HoKo->family[n_cluster]<<endl;
        return 0;
    }
    
    //This variable is used to assing node numbers and after assigning all node numbers
    //it will contain the number of nodes in the network
    //It is initialized with the first available node, which is equal to reserved_nodes
    long int global_nodes = (long int)reserved_nodes;
    
    //Construct the LM matrices depending on the particles available on each cluster
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) {
        
        //There are CNT clusters, so construct the LM matrix for CNTs
        //hout<<"LM_matrix_for_cnts"<<endl;
        if (!LM_matrix_for_cnts(n_cluster, HoKo, Cutwins, global_nodes)) {
            hout<<"Error in Compute_voltage_field when calling LM_matrix_for_cnts"<<endl;
            return 0;
        }
    }
    if (HoKo->clusters_gnp.size() && HoKo->clusters_gnp[n_cluster].size()) {
        
        //There are GNP clusters, so construct the LM matrix for GNPs
        //hout<<"LM_matrix_for_gnps"<<endl;
        if (!LM_matrix_for_gnps(n_cluster, HoKo, Cutwins, structure_gnp, global_nodes)) {
            hout<<"Error in Compute_voltage_field when calling LM_matrix_for_gnps"<<endl;
            return 0;
        }
    }

    //hout << "global_nodes="<<global_nodes<<endl;
    //Variables for using the SSS for the sparse matrix
    vector<long int> col_ind, row_ptr;
    vector<double> values, diagonal(global_nodes, 0);
    //P is the search direction
    //R is the residual vector
    vector<double> P(global_nodes-reserved_nodes,0), R(global_nodes-reserved_nodes,0);
    //Prescribed voltages applied to the sample
    vector<double> VEF;
    
    //With the LM matrix, now fill the sparse stiffness matrix
    //Fill_sparse_stiffness_matrix (use R_flag)
    hout<<"Fill_sparse_stiffness_matrix"<<endl;
    if (!Fill_sparse_stiffness_matrix(R_flag, global_nodes, reserved_nodes, cutoffs.van_der_Waals_dist, n_cluster, electric_param, HoKo, points_cnt, radii, structure_gnp, points_gnp, gnps, col_ind, row_ptr, values, diagonal, P, R, VEF)) {
        hout<<"Error in Compute_voltage_field when calling Fill_sparse_stiffness_matrix"<<endl;
        return 0;
    }
    
    hout << "Solve_DEA_equations_CG_SSS"<<endl;
    //This is where the actual direct electrifying algorithm (DEA) takes place
    if (!Solve_DEA_equations_CG_SSS(global_nodes, reserved_nodes, col_ind, row_ptr, values, diagonal, P, R, VEF)) {
        hout<<"Error in Compute_voltage_field when calling Solve_DEA_equations_CG_SSS"<<endl;
        return 0;
    }
    
    return 1;
}
//This function gets the first available node depending on the percolated family
int Direct_Electrifying::Get_global_nodes(const int &family)
{
    if (family == 6) {
        //If family is 6, then a cluster percolates in the three directions.
        //Hence we need 6 voltages: 0 to 5.
        //6 is the first available node
        return 6;
    } else if ( 3 <= family && family <= 5 ){
        //If family is 3, 4 or 5, then a cluster percolates in two directions.
        //Hence we need 4 voltages: 0 to 3.
        //4 is the first available node
        return 4;
    } else if ( 0 <= family && family <= 2 ) {
        //If a family is 0, 1 or 2, then a cluster percolates in one direction.
        //Hence we need 2 voltages: 0 to 1.
        //2 is the first available node
        return 2;
    }
    else {
        //There was an error, this clause should not be reached
        return -1;
    }
}
//This function fills the LM matrix for CNTs
int Direct_Electrifying::LM_matrix_for_cnts(const int &n_cluster, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, long int &global_nodes)
{
    //Map all CNT points at a boundary with prescribed boundary conditions
    //First use a temporary variable to save time when mappong the first and
    //last points of a CNT
    map<long int, long int> LMM_cnts_boundary;
    //hout<<"Map_points_at_boundaries"<<endl;
    if (!Map_points_at_boundaries(HoKo->family[n_cluster], Cutwins->boundary_cnt_pts, LMM_cnts_boundary)) {
        hout<<"Error in LM_matrix_for_cnts when calling Map_points_at_boundaries (1)"<<endl;
        return 0;
    }
    
    //Iterate over all CNTs in the cluster
    for (int i = 0; i < (int)HoKo->clusters_cnt[n_cluster].size(); i++) {
        
        //Get the current CNT
        int CNT = HoKo->clusters_cnt[n_cluster][i];
        //hout<<"HoKo->elements_cnt[CNT="<<CNT<<"].size="<<HoKo->elements_cnt[CNT].size()<<endl;
        
        //Iterator for the elements in set
        set<long int>::iterator j = HoKo->elements_cnt[CNT].begin();
        
        //Iterator to find an element in LMM_cnts_boundary
        map<long int, long int>::iterator it;
        
        //Chek if the first point in the CNT element is at a boundary
        //hout<<"LMM_cnts_boundary.find("<<*j<<") == LMM_cnts_boundary.end()"<<endl;
        if (LMM_cnts_boundary.find(*j) == LMM_cnts_boundary.end()) {
            
            //Initial point of CNT is not at a boundary, so map it to a new node
            LMM_cnts[*j] = global_nodes;
            
            //Update the next available node
            global_nodes++;
        }
        
        //Iterator at the last element of the set
        set<long int>::iterator j_end = HoKo->elements_cnt[CNT].end();
        //end() points after the last element, so decrease the iterator to point after
        //the penultimate element and thus to the last element
        j_end--;
        //hout<<"j_end="<<*j_end<<" rbegin="<<*(HoKo->elements_cnt[CNT].rbegin())<<endl;
        
        //Map all points in the elements vector
        //hout<<"Map all points in the elements vector"<<endl;
        for (j++; j != j_end; j++) {
            
            //Get the current point number
            long int P = *j;
            //hout<<"P="<<P<<" global_nodes="<<global_nodes<<endl;
            
            //Map the point to a new node
            LMM_cnts[P] = global_nodes;
            
            //Update the next available node
            global_nodes++;
        }
        
        //Chek if the last point in the CNT element is at a boundary
        //hout<<"LMM_cnts_boundary.find("<<*j_end<<") == LMM_cnts_boundary.end()"<<endl;
        if (LMM_cnts_boundary.find(*j_end) == LMM_cnts_boundary.end()) {
            
            //Last point of CNT is not at a boundary, so map it to a new node
            LMM_cnts[*j_end] = global_nodes;
            
            //Update the next available node
            global_nodes++;
        }
    }
    
    //Add the mappings of all CNT points at a boundary with prescribed boundary conditions
    //into the class variable
    //hout<<"Map_points_at_boundaries"<<endl;
    if (!Map_points_at_boundaries(HoKo->family[n_cluster], Cutwins->boundary_cnt_pts, LMM_cnts)) {
        hout<<"Error in LM_matrix_for_cnts when calling Map_points_at_boundaries (2)"<<endl;
        return 0;
    }
    
    return 1;
}
//
int Direct_Electrifying::Map_points_at_boundaries(const int &family, const vector<vector<long int> > &boundary_pts, map<long int, long int> &LMM)
{
    //Get the vector of boundaries
    vector<int> boundaries;
    if (!Get_vector_of_boundaries(family, boundaries)) {
        hout<<"Error in Map_cnt_points_at_boundaries when calling Get_vector_of_boundaries"<<endl;
        return 0;
    }
    
    //Iterate over the vector of boundaries to fill the map LMM_cnts
    //The index of boundaries is the node number with presecribed boundary conditions
    for (long int n = 0; n < (int)boundaries.size(); n++) {
        
        //Get the boundary number
        int b = boundaries[n];
        
        //Iterate over the points at boundary b
        for (int i = 0; i < (int)boundary_pts[b].size(); i++) {
            
            //Get the current point number
            long int P = boundary_pts[b][i];
            
            //Add the mapping from P to node n
            LMM[P] = n;
        }
    }
    
    return 1;
}
//This function generates a vector of boundaries to be scanned for the local mapping of
//reserved nodes
int Direct_Electrifying::Get_vector_of_boundaries(const int &family, vector<int> &boundaries)
{
    //Add boundary numbers to the vector boundaries based on the family
    switch (family) {
        //Family X
        case 0:
            boundaries.push_back(2);
            boundaries.push_back(4);
            break;
        //Family Y
        case 1:
            boundaries.push_back(3);
            boundaries.push_back(5);
            break;
        //Family Z
        case 2:
            boundaries.push_back(0);
            boundaries.push_back(1);
            break;
        //Family XY
        case 3:
            boundaries.push_back(2);
            boundaries.push_back(4);
            boundaries.push_back(3);
            boundaries.push_back(5);
            break;
        //Family XZ
        case 4:
            boundaries.push_back(2);
            boundaries.push_back(4);
            boundaries.push_back(0);
            boundaries.push_back(1);
            break;
        //Family YZ
        case 5:
            boundaries.push_back(3);
            boundaries.push_back(5);
            boundaries.push_back(0);
            boundaries.push_back(1);
            break;
        //Family XYZ
        case 6:
            boundaries.push_back(2);
            boundaries.push_back(4);
            boundaries.push_back(3);
            boundaries.push_back(5);
            boundaries.push_back(0);
            boundaries.push_back(1);
            break;
            
        default:
            hout<<"Error in Get_vector_of_boundaries, invalid family: "<<family<<endl;
            return 0;
    }
    
    return 1;
}
//This function fills the LM matrix for GNPs
int Direct_Electrifying::LM_matrix_for_gnps(const int &n_cluster, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure_gnp, long int &global_nodes)
{
    
    //Add the mappings of all GNP points at a boundary with prescribed boundary conditions
    //into the class variable
    if (!Map_points_at_boundaries(HoKo->family[n_cluster], Cutwins->boundary_gnp_pts, LMM_gnps)) {
        hout<<"Error in LM_matrix_for_gnps when calling Map_points_at_boundaries"<<endl;
        return 0;
    }
    
    //Iterate over all GNPs in the cluster
    for (int i = 0; i < (int)HoKo->clusters_cnt[n_cluster].size(); i++) {
        
        //Get current GNP
        int GNP = HoKo->clusters_cnt[n_cluster][i];
        
        //Iterate over all points in GNP
        for (int j = 0; j < (int)structure_gnp[GNP].size(); j++) {
            
            //Get current point number
            long int P = structure_gnp[GNP][j];
            
            //Check it is not a boundary node
            if (P >= Cutwins->n_gnp_pts_b) {
                
                //P is not a boundary point
                
                //Add the mapping of the current GNP point to a node number
                LMM_gnps[P] = global_nodes;
                
                //Update the next available node
                global_nodes++;
            }
        }
    }
    
    return 1;
}
//This function creates the sparse stifness matrix that will be used to solve the sytem of equations
//The sparse version is more efficient computationally speaking
int Direct_Electrifying::Fill_sparse_stiffness_matrix(const int &R_flag, const long int &nodes, const long int &reserved_nodes, const double &d_vdw, const int &n_cluster, const Electric_para &electric_param, Hoshen_Kopelman *HoKo, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, vector<double> &P, vector<double> &R, vector<double> &VEF)
{
    //------------------------------------------------------------------------
    //Initially, generate a 2D vector that will be used to generate the 1D vectors for SSS
    //col_values[i] has the map col_j->value_ij, where i is the row number
    vector<map<long int, double>> col_values(nodes);
    
    //------------------------------------------------------------------------
    //Add contributions from particles
    
    //Fill the 2D matrices with the contributions of the CNTs when there are CNT clusters
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) {
        
        //Add contributions from CNT resistors
        //hout << "Fill_2d_matrices_cnts"<<endl;
        if (!Fill_2d_matrices_cnts(R_flag, n_cluster, electric_param, HoKo, points_cnt, radii, col_values, diagonal)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_cnts" << endl;
            return 0;
        }
        
        //Add contributions from CNT-CNT junctions if any
        if (HoKo->cluster_cnt_junctions.size() && HoKo->cluster_cnt_junctions[n_cluster].size()) {
            //hout << "Fill_2d_matrices_cnt_junctions"<<endl;
            if (!Fill_2d_matrices_cnt_junctions(R_flag, d_vdw, electric_param, HoKo->cluster_cnt_junctions[n_cluster], HoKo->junctions_cnt, points_cnt, radii, LMM_cnts, col_values, diagonal)) {
                hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_cnt_junctions" << endl;
                return 0;
            }
        }
    }
    
    //Set used to determine the nodes from CNT point in mixed junctions
    //map<long int, long int> points_cnt_rad;
    
    //Check if there are mixed junctions for n_cluster
    if (HoKo->cluster_mix_junctions.size() && HoKo->cluster_mix_junctions[n_cluster].size()) {
        
        //Add contributions from mixed junctions
        //hout << "Fill_2d_matrices_mixed_junctions"<<endl;
        if (!Fill_2d_matrices_mixed_junctions(R_flag, d_vdw, electric_param, HoKo->cluster_mix_junctions[n_cluster], HoKo->junctions_mixed, points_cnt, radii, points_gnp, gnps, LMM_cnts, LMM_gnps, col_values, diagonal, points_cnt_rad)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_mixed_junctions" << endl;
            return 0;
        }
    }
    
    //Fill the 2D matrices with the contributions of the GNPs when there are GNP clusters
    if (HoKo->clusters_gnp.size() && HoKo->clusters_gnp[n_cluster].size()) {
        
        //Add contributions from GNP resistors
        //hout<<"Fill_2d_matrices_gnp"<<endl;
        if (!Fill_2d_matrices_gnp(R_flag, electric_param, HoKo->clusters_gnp[n_cluster], points_gnp, structure_gnp, gnps, LMM_gnps, points_cnt_rad, col_values, diagonal)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_gnp" << endl;
            return 0;
        }
        
        //Add contributions from GNP-GNP junctions, if any
        //hout<<"Fill_2d_matrices_cnts"<<endl;
        if (HoKo->cluster_gnp_junctions.size() && HoKo->cluster_gnp_junctions[n_cluster].size()) {
            if (!Fill_2d_matrices_gnp_junctions(R_flag, d_vdw, electric_param, HoKo->cluster_gnp_junctions[n_cluster], HoKo->junctions_gnp, points_gnp, gnps, LMM_gnps, col_values, diagonal)) {
                hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_gnp_junctions" << endl;
                return 0;
            }
        }
    }
    
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------
    //Convert from 2D vectors to 1D vectors
    
    //Initialize row_ptr
    //The first element of row_ptr is zero
    row_ptr.push_back(0);
    vector<vector<double> > KEFT(nodes-reserved_nodes, vector<double> (nodes,0));
    
    //hout<<"From_2d_to_1d_vectors"<<endl;
    if (!From_2d_to_1d_vectors(reserved_nodes, col_values, KEFT, col_ind, row_ptr, values, diagonal)) {
        hout<<"Error in Fill_sparse_stiffness_matrix when calling From_2d_to_1d_vectors"<<endl;
        return 0;
    }
    
    //=========================================
    //Set up variables for the Conjugate Gradient Algorithm:
    //Search direction (P) and residual vector (R)
    //hout<<"Set_up_residual_and_search_direction"<<endl;
    if (!Set_up_residual_and_search_direction(R_flag, nodes, reserved_nodes, electric_param, KEFT, P, R, VEF)) {
        hout<<"Error in Fill_sparse_stiffness_matrix when calling Set_up_residual_and_search_direction"<<endl;
        return 0;
    }
    
    //Print some vectors
    Printer Pr;
    Pr.Print_1d_vec(diagonal, "diagonal.txt");
    Pr.Print_1d_vec(row_ptr, "row_ptr.txt");
    Pr.Print_1d_vec(col_ind, "col_ind.txt");
    Pr.Print_1d_vec(values, "values.txt");
    
    return 1;
}
//This function adds the contributions of the CNT and junction resistors to the stiffness matrix
int Direct_Electrifying::Fill_2d_matrices_cnts(const int &R_flag, const int &n_cluster, const Electric_para &electric_param,  Hoshen_Kopelman *HoKo, const vector<Point_3D> &points_cnt, const vector<double> &radii, vector<map<long int, double>> &col_values, vector<double> &diagonal)
{
    //Scan every CNT in the cluster
    for (int i = 0; i < (long int)HoKo->clusters_cnt[n_cluster].size(); i++) {
        
        //Current CNT in cluster
        int CNT = HoKo->clusters_cnt[n_cluster][i];
        
        //Iterator for looping over the elements in set
        set<long int>::const_iterator it = HoKo->elements_cnt[CNT].begin();
        
        //Get the first point of the element, i.e., the first point of the CNT
        long int P1 = *it;
        
        //Iterate over the points in the elements vector of CNT
        for (it++; it != HoKo->elements_cnt[CNT].end(); it++) {
            
            //Get the second point of the element
            long int P2 = *it;
            
            //Calculate resistance depending of the R_flag
            double Re;
            if (!Calculate_resistance_cnt(R_flag, points_cnt, P1, P2, radii[CNT], electric_param.resistivity_CNT, Re)) {
                hout<<"Error in Fill_2d_matrices_cnts when calling Calculate_resistance_cnt"<<endl;
                return 0;
            }
            
            //Get the nodes of the points
            long int node1 = LMM_cnts[P1];
            long int node2 = LMM_cnts[P2];
            
            //Add the resistance value to the corresponding 2D vectors
            if(!Add_new_elements_to_2d_sparse_matrix(node1, node2, Re, col_values, diagonal)) {
                hout<<"Error in Fill_2d_matrices_cnts when calling Add_elements_to_2d_sparse_matrix"<<endl;
                return 0;
            }
            
            //Set the second point as the first point for the next iteration
            P1 = P2;
        }
    }
    
    return 1;
}
//This function calculates the length of a CNT segment that corresponds to one element (resistor)
//Using the resisitivity as input parameter the resistance of the CNT segment is calculated
int Direct_Electrifying::Calculate_resistance_cnt(const int &R_flag, const vector<Point_3D> &points_cnt, const long int &P1, const long int &P2, const double &radius, const double &resistivity, double &Re)
{
    //Check the resistance flag
    if (R_flag == 1) {
        
        //Calculate the actual resistance
        
        //Variable to store the CNT length
        double length = 0;
        //Calculate the length of each CNT segment
        for (long int i = P1; i < P2; i++) {
            //Calculate the distance from point i to i+1 and add it to the total length
            length = length + points_cnt[i].distance_to(points_cnt[i+1]);
        }
        
        //Calculate resistance as R = rho*l/A; A = PI*r*r
        Re = resistivity*length/(PI*radius*radius);
    }
    else if (R_flag == 0) {
        
        //Use unit value
        Re = 1.0;
    }
    else {
        hout << "Error in Fill_2d_matrices_cnts. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
        return 0;
    }
    
    return 1;
}
//
int Direct_Electrifying::Add_new_elements_to_2d_sparse_matrix(const long int &node1, const long int &node2, const double &Re, vector<map<long int, double> > &col_values, vector<double> &diagonal)
{
    double Re_inv = 1/Re;
    //Add the diagonal elements of the stiffness matrix
    diagonal[node1] += Re_inv;
    diagonal[node2] += Re_inv;
    
    //Add the off diagonal elements of the stiffness matrix
    if (node1 > node2) {
        
        //node2 is the column index and node1 is the row index
        //In row node1, add the map node2-> -Re_inv
        col_values[node1][node2] = -Re_inv;
    } else {
        
        //node1 is the column index and node2 is the row index
        //In row node2, add the map node1-> -Re_inv
        col_values[node2][node1] = -Re_inv;
    }
    //hout << "Added ";
    return 1;
}
//This function checks if an element was already added to the 2D sparse stiffness matrix
//If already in the matrix, the contribution of the new element is added to the existing value
int Direct_Electrifying::Add_to_existing_elements_in_2d_sparse_matrix(const long int &node1, const long int &node2, const double &Re, vector<map<long int, double> > &col_values, vector<double> &diagonal)
{
    double Re_inv = 1/Re;
    //Add the diagonal elements of the stiffness matrix
    diagonal[node1] += Re_inv;
    diagonal[node2] += Re_inv;
    
    //Add the off diagonal elements of the stiffness matrix
    if (node1 > node2) {
        
        //node2 is the column index and node1 is the row index
        //In row node1, add the map node2-> -Re_inv
        if (col_values[node1].find(node2) == col_values[node1].end()) {
            
            //node2 is not a key in the map col_values[node1], so just add a new map
            col_values[node1][node2] = -Re_inv;
        }
        else {
            
            //node2 is already a key in the map col_values[node1], so add to existing value
            col_values[node1][node2] = col_values[node1][node2] -Re_inv;
        }
        
    } else {
        
        //node1 is the column index and node2 is the row index
        //In row node2, add the map node1-> -Re_inv
        if (col_values[node2].find(node1) == col_values[node2].end()) {
            
            //node1 is not a key in the map col_values[node2], so just add a new map
            col_values[node2][node1] = -Re_inv;
        }
        else {
            
            //node2 is already a key in the map col_values[node1], so add to existing value
            col_values[node2][node1] = col_values[node2][node1] -Re_inv;
        }
    }
    //hout << "Added ";
    return 1;
}
int Direct_Electrifying::Fill_2d_matrices_cnt_junctions(const int &R_flag, const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_cnt_junctions_i, const vector<Junction> &junctions_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const map<long int, long int> &LMM_cnts, vector<map<long int, double> > &col_values, vector<double> &diagonal)
{
    //Iterate over all junctions in the cluster
    //hout<<"cluster_cnt_junctions_i.size()="<<cluster_cnt_junctions_i.size()<<endl;
    for (int i = 0; i < (int)cluster_cnt_junctions_i.size(); i++) {
        
        //Get current junction index
        int idx = cluster_cnt_junctions_i[i];
        
        //Get the point numbers of the mixed juntion
        long int P1 = junctions_cnt[idx].P1;
        long int P2 = junctions_cnt[idx].P2;
        //hout<<"P1="<<P1<<" P2="<<P2<<endl;
        
        //Get the particle numbers
        int CNT1 = junctions_cnt[idx].N1;
        //hout<<"CNT1="<<CNT1<<endl;
        int CNT2 = junctions_cnt[idx].N2;
        //hout<<"CNT2="<<CNT2<<endl;
        
        //Calculate the junction resistance
        double Re = 1.0;
        if (R_flag == 1) {
            if (!Calculate_junction_resistance(junctions_cnt[idx], d_vdw, radii[CNT1], points_cnt[P1], radii[CNT2], points_cnt[P2], electric_param, Re)) {
                hout<<"Error in Fill_2d_matrices_cnt_junctions when calling Calculate_junction_resistance"<<endl;
                return 0;
            }
        }
        else if (R_flag != 0) {
            hout << "Error in Fill_2d_matrices_cnt_junctions. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
            return 0;
        }
        
        //Get node numbers
        long int node1 = LMM_cnts.at(P1);
        long int node2 = LMM_cnts.at(P2);
        //hout<<"node1="<<node1<<" node2="<<node2<<endl;
        
        //Add junction resistance to sparse stiffness matrix
        //hout<<"Add_new_elements_to_2d_sparse_matrix i="<<i<<endl;
        if (!Add_new_elements_to_2d_sparse_matrix(node1, node2, Re, col_values, diagonal)) {
            hout<<"Error in Fill_2d_matrices_cnt_junctions when calling Add_elements_to_2d_sparse_matrix"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function calculates the junction resistance between either two CNTs, two GNPs, or
//a CNT and a GNP
int Direct_Electrifying::Calculate_junction_resistance(const Junction &j, const double &d_vdw, const double &rad1, const Point_3D &P1, const double &rad2, const Point_3D &P2, const struct Electric_para &electric_param, double &Re)
{
    //Check if a constant resistance is used
    if (electric_param.junction_type == "constant") {
        Re = electric_param.junction_resistance;
        return 1;
    }
    else if (electric_param.junction_type == "exponential") {
        
        //Calcualte the distance between the to points
        double separation = P1.distance_to(P2);
        
        //If the first particle of the junction is a CNT, susbtract the radius
        if (j.type1 == "CNT") {
            separation = separation - rad1;
        }
        
        //If the second particle of the junction is a CNT, susbtract the radius
        if (j.type2 == "CNT") {
            separation = separation - rad2;
        }
        
        //Check if the tunnel distance is below the van der Waals distance
        //This happens when the penetrating model is used
        if (separation - d_vdw < Zero) {
            
            //If the distance between the points is below the van der Waals distance, then set the separation equal to the van der Waals distance
            separation = d_vdw;
        }
        
        //Calculate the area of the junction assuming a circle using the smalles radius
        //which will the one limiting the flow of electrons
        double A = (rad1 - rad2 < Zero)? PI*rad1*rad1: PI*rad2*rad2;
        
        //Calculate quantity associted with a squared root
        double sqrt_tmp = sqrt(2*electric_param.e_mass*electric_param.lambda_barrier*electric_param.e_charge);
        
        //Calculate the exponential term
        double exp_tmp = exp(4000*PI*separation*sqrt_tmp/electric_param.h_plank);
        
        //Calculate term that multiplies the exponential
        double denominator_tmp = A*electric_param.e_charge*electric_param.e_charge*sqrt_tmp;
        double mult_tmp = 10*electric_param.h_plank*electric_param.h_plank*separation/denominator_tmp;
        
        //Calculate tunnel resistance
        Re = mult_tmp*exp_tmp;
    }
    else {
        hout<<"Error in Calculate_junction_resistance: invalid junction type. Input was: "<<electric_param.junction_type<<endl;
        return 0;
    }
    
    return 1;
}
//This function adds the contributions of the junctions between a CNT and a GNP
int Direct_Electrifying::Fill_2d_matrices_mixed_junctions(const int &R_flag, const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_mix_junctions_i, const vector<Junction> &junctions_mixed, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const map<long int, long int> &LMM_cnts, const map<long int, long int> &LMM_gnps, vector<map<long int, double> > &col_values, vector<double> &diagonal, map<long int, double> &points_cnt_rad)
{
    //Iterate over all junctions in the cluster
    for (int i = 0; i < (int)cluster_mix_junctions_i.size(); i++) {
        
        //Get current junction index
        int idx = cluster_mix_junctions_i[i];
        
        //Get the point numbers of the mixed juntion
        long int Pcnt = junctions_mixed[idx].P1;
        long int Pgnp = junctions_mixed[idx].P2;
        
        //Get the particle numbers
        int cnt_n = junctions_mixed[idx].N1;
        int gnp_n = junctions_mixed[idx].N2;
        
        //Calculate the junction resistance
        double Re = 1.0;
        if (R_flag == 1) {
            if (!Calculate_junction_resistance(junctions_mixed[idx], d_vdw, radii[cnt_n], points_cnt[Pcnt], gnps[gnp_n].t/2, points_gnp[Pgnp], electric_param, Re)) {
                hout<<"Error in Fill_2d_matrices_mixed_junctions when calling Calculate_junction_resistance"<<endl;
                return 0;
            }
        }
        else if (R_flag != 0) {
            hout << "Error in Fill_2d_matrices_mixed_junctions. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
            return 0;
        }
        
        //Get node numbers
        long int node1 = LMM_cnts.at(Pcnt);
        long int node2 = LMM_gnps.at(Pgnp);
        
        //Add junction resistance to sparse stiffness matrix
        if (!Add_new_elements_to_2d_sparse_matrix(node1, node2, Re, col_values, diagonal)) {
            hout<<"Error in Fill_2d_matrices_mixed_junctions when calling Add_elements_to_2d_sparse_matrix"<<endl;
            return 0;
        }
        
        //Add the radius of the CNT that has point P into the set of nodes with CNT radius
        points_cnt_rad[Pgnp] = radii[cnt_n];
    }
    
    return 1;
}
//Add contributions from the resistors formed by the triangulation on the GNP
int Direct_Electrifying::Fill_2d_matrices_gnp(const int &R_flag, const Electric_para &electric_param, const vector<int> &cluster_gnp, const vector<Point_3D> &points_gnp, const vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, const map<long int, long int> &LMM_gnps, map<long int, double> &points_cnt_rad, vector<map<long int, double> > &col_values, vector<double> &diagonal)
{
    //Triangulation object
    Triangulation dt;
    
    //Scan every GNP, perform the triangulations and add the elements to the stiffness matrix
    for (long int i = 0; i < (long int)cluster_gnp.size(); i++) {
        
        //current hybrid particle
        int gnp_i = cluster_gnp[i];
        
        //Check if the triangulation needs to be performed
        //Triangulation is performed when calculating the backbone (R_flag == 0)
        //or when calculating the electrical resistance and there is already a triangulation
        if (!R_flag || gnps[gnp_i].triangulation.empty()) {
            
            //Perform triangulation
            if (!dt.Generate_3d_trangulation(points_gnp, structure_gnp[gnp_i], gnps[gnp_i])) {
                hout << "Error in Fill_2d_matrices_gnp when calling delaunay->Generate_3d_trangulation" << endl;
                return 0;
            }
            //hout << "Triangulation " << i << ", " << structure_gnp[GNP].size() << " points in GNP " << GNP<<", " << gnps[GNP].triangulation.size() << " edges"<< endl;
        }
        
        //Add elements from the triangulation
        //Scan triangulation edges in reverse order since some edges might need to be deleted
        for (int j = (int)gnps[gnp_i].triangulation.size() - 1; j >= 0 ; j--) {
            
            //Get the two vertices of the triangulation
            long int v1 = gnps[gnp_i].triangulation[j].v1;
            long int v2 = gnps[gnp_i].triangulation[j].v2;
            
            //Get the node numbers
            long int node1 = LMM_gnps.at(v1);
            long int node2 = LMM_gnps.at(v2);
            
            //Initialize resistor with unit resistance
            double Re = 1.0;
            
            //Check if the actual resistance needs to be calculated
            if (R_flag == 1) {
                
                //Radii for the junction ressitance
                double rad1 = 0, rad2 = 0;
                
                //Check if v1 has a junction with a CNT or is a CNT seed
                if (points_cnt_rad.find(v1) == points_cnt_rad.end()) {
                    //v1 comes from a GNP-GNP junction so use half the thickness as rad1
                    rad1 = gnps[gnp_i].t/2;
                }
                else {
                    //v1 comes from a CNT-GNP junction or it is a CNT seed, so use CNT radius
                    //Get radius
                    rad1 = points_cnt_rad[v1];
                }
                
                //Check if v2 has a junction with a CNT or is a CNT seed
                if (points_cnt_rad.find(v2) == points_cnt_rad.end()) {
                    //v2 comes from a GNP-GNP junction so use half the thickness as rad2
                    rad2 = gnps[gnp_i].t/2;
                }
                else {
                    //v1 comes from a CNT-GNP junction or it is a CNT seed, so use CNT radius
                    //Get radius
                    rad2 = points_cnt_rad[v2];
                }
                
                //Calculate the triangulation resistor, i.e., the resistance of
                //the "conduction band" in the GNP
                if (!Calculate_resistance_gnp(points_gnp[v1], points_gnp[v2], rad1, rad2, electric_param, Re)) {
                    hout << "Error in Fill_2d_matrices_gnp when calling Calculate_resistance_gnp" << endl;
                    return 0;
                }
                
            }
            else if (R_flag != 0) {
                hout << "Error in Fill_2d_matrices_gnp. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
                return 0;
            }
            
            if (!Add_new_elements_to_2d_sparse_matrix(node1, node2, Re, col_values, diagonal)) {
                hout << "Error in Fill_2d_matrices_gnp when calling Add_elements_to_2d_sparse_matrix" << endl;
                return 0;
            }
            //hout << "flag1="<<hybrid_particles[hyb].triangulation_flags[j][0]<<" P1="<<P1<<" node1="<<node1<<" flag2="<<hybrid_particles[hyb].triangulation_flags[j][1]<<" P2="<<P2<<" node2="<<node2<<" Re="<<Re<<endl;
        }
    }
    //hout << "Triangulation done" << endl;
    
    //
    //NOTE:
    //When adding hybrid particles, the stiffness matrix need to be filled differently
    //Thus, a new function is needed for hybrid particles due to the special cases that
    //arise from having CNT points attached to GNPs
    //This if-statement is left commented, but when adding hybrid particles into the code
    //it needs to be used
    /*
    if (particle_type == "Hybrid_particles") {
        for (long int i = 0; i < (long int)cluster_gnp.size(); i++) {
            //New function that consideres the cases below:
            //CASE 1: Check if the edge has to be ignored
            //Sometimes a CNT seed is also a boundary point, in that case the triangulation
            //will result in having repeated elements in the stiffness matrix
            //This since two boundary nodes have the same node number, thus this adds a
            //resistor connected to the same node on both ends
            //Because of this, tunneling is ignored on a bonudary point where voltage is applied
            //Thus ignore elements that have two points that are boundary nodes
            if (node1 == node2 && node1 < reserved_nodes) {
                //The condition "node2 < reserved_nodes" is not needed
                //It is redundant because of the first condition, if node1 and node 2 are equal
                //and node1 is less than reserved_nodes, then for sure node2 is also
                //less than reserved_nodes
                
                //If the edge needs to be ignored in the construction of the stiffness matrix
                //Delete the edge
                gnps[GNP].triangulation.erase(gnps[GNP].triangulation.begin()+j);
            }
            else {
                
                //Whatever the case is, the resistance is needed
                //So first calculate the resistance of the edge depending on the value of R_flag
     
                //Check if the resistor radius should be a CNT radius (hybrid or mixed
                //particles) or the GNP thicknes (only GNPs)
                //This is done by choosing the minimum point
                
                //CASE 2:
                //Sometimes a triangulation edge has a CNT seed connected to a boundary GNP point without any other node in between
                //In this case, the CNT resistor and the triangulation edge connect the same nodes so this results
                //in repeated elements in the stiffness matrix
                //Thus, check if this is happening
                if (CASE2) {
                    //Deal with case 2
                }
                else {
                    //CASE 3:
                    //Add resistor following the standard procedure (code within the loop below)
                }
            }
        }
    }
    else {
        //Add all resistors following the standadard procedure (the code below)
    }*/
    
    return 1;
}
//This function calculates the resistance that comes from the triangulation on the GNPs
//It uses the resistivities along the surface and along the thickness so the resistance
//can be calculated along any direction the triangulation edge may have
//Flag:
//0: CNT-CNT edge
//1: GNP-GNP edge
//2: CNT-GNP edge
int Direct_Electrifying::Calculate_resistance_gnp(const Point_3D &P1, const Point_3D &P2, const double &rad1, const double &rad2, const struct Electric_para &electric_param, double &Re)
{
    //Calculate the distance between the points in contact
    double L = P1.distance_to(P2);
    
    //Unit vector in the direction from P1 to P2
    Point_3D u_direction = (P2 - P1)/L;
    
    //Dot product of resistivity tensor with the u_direction vector
    //i.e., rho.dot(u_direction)
    Point_3D resistivity(u_direction.x*electric_param.resistivity_GNP_surf, u_direction.y*electric_param.resistivity_GNP_surf, u_direction.z*electric_param.resistivity_GNP_t);
    
    //Get the dot product of the resistivity vector with the direction vector
    //i.e., u_direction.dot(rho.dot(u_direction))
    double rho = resistivity.dot(u_direction);
    
    //Calcualte resistance as:
    //Re = rho*L/(PI*rad1*rad2)
    Re = rho*L/(PI*rad1*rad2);
    
    return 1;
}
//This function adds GNP-GNP junctions to the 2D sparse stiffness matrix
int Direct_Electrifying::Fill_2d_matrices_gnp_junctions(const int &R_flag, const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_gnp_junctions_i, const vector<Junction> &junctions_gnp, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const map<long int, long int> &LMM_gnps, vector<map<long int, double> > &col_values, vector<double> &diagonal)
{
    //Iterate over all junctions in the cluster
    for (int i = 0; i < (int)cluster_gnp_junctions_i.size(); i++) {
        
        //Get current junction index
        int idx = cluster_gnp_junctions_i[i];
        
        //Get the point numbers of the mixed juntion
        long int P1 = junctions_gnp[idx].P1;
        long int P2 = junctions_gnp[idx].P2;
        
        //Get the particle numbers
        int GNP1 = junctions_gnp[idx].N1;
        int GNP2 = junctions_gnp[idx].N2;
        
        //Calculate the junction resistance
        double Re = 1.0;
        if (R_flag == 1) {
            if (!Calculate_junction_resistance(junctions_gnp[idx], d_vdw, gnps[GNP1].t/2, points_gnp[P1], gnps[GNP2].t/2, points_gnp[P2], electric_param, Re)) {
                hout<<"Error in Fill_2d_matrices_gnp_junctions when calling Calculate_junction_resistance"<<endl;
                return 0;
            }
        }
        else if (R_flag != 0) {
            hout << "Error in Fill_2d_matrices_gnp_junctions. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
            return 0;
        }
        
        //Get node numbers
        long int node1 = LMM_gnps.at(P1);
        long int node2 = LMM_gnps.at(P2);
        
        //Add junction resistance to sparse stiffness matrix
        if (!Add_new_elements_to_2d_sparse_matrix(node1, node2, Re, col_values, diagonal)) {
            hout<<"Error in Fill_2d_matrices_mixed_junctions when calling Add_elements_to_2d_sparse_matrix"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function transforms the 2D vectors that contain the stiffness matrix into 1D vectors
//so they can be in the SSS format and make the matrix-vector multiplications faster
//For reference, the 2D stiffness matrix has this form:
//
//   | KE    KEF |
//   | KEFT  KF  |
//
//For the CG algorithm we need to extract KEFT and KF.
//We actually already have the lower left corner of the stiffness matrix
//
int Direct_Electrifying::From_2d_to_1d_vectors(const long int &reserved_nodes, const vector<map<long int, double> > &col_values, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal)
{
    //hout << "Fill 1D vectors" <<endl;
    //Iterate over the rows in the 2D vector of the stiffness matrix (col_values)
    //Skip the top rows of the matrix that correspond to the reserved nodes
    for (long int i = reserved_nodes; i < (long int)col_values.size(); i++) {
        
        //hout << "for(i)=" <<i << " col_values[i].size()="<<col_values[i].size()<< endl;
        //Iterate over the columns in row i
        for (map<long int, double>::const_iterator it = col_values[i].begin(); it != col_values[i].end(); it++) {
            
            //Get the column number
            long int col = it->first;
            
            //Get the value in row i, column col
            double Kij = it->second;
            //hout<<"\tcol="<<col<<" Kij="<<Kij<<endl;
            
            //Check if the column corresponds to a reserved node
            if (col >= reserved_nodes) {
                
                //When the colum index is greater or equal to reserved_nodes,
                //then it is the lower-righ of the full stiffness matrix
                values.push_back(Kij);
                
                //The column numbers that I need are the numbers in the full matrix-reserved_nodes
                col_ind.push_back(col - reserved_nodes);
            }
            else {
                
                //When the column index is less than reserved_nodes,
                //I need to save the value of the resistance on the vector KEFT
                //so that it is used for the CG algorithm.
                //The column index is (of course) the column index of KEFT as
                //in this else-statement we have that 0<=col<reserved_nodes
                //The row is the index iterator i-reserved_nodes
                KEFT[i-reserved_nodes][col] = Kij;
            }
        }
        
        //hout << "row_ptr" << endl;
        //Add the number of non-zero entries in the lower tirangle of the stiffness matrix
        //that have been added to the vector values
        row_ptr.push_back(values.size());
    }
    
    //Remove the elements of the diagonal that are not used
    //hout << "Remove the elements of the diagonal that are not used"<<endl;
    diagonal.erase(diagonal.begin(), diagonal.begin()+reserved_nodes);
    
    return 1;
}
//This functio sets up the search direction (P) and residual vector (R) for the CG, and also
//the vector of prescribed voltages VEF
int Direct_Electrifying::Set_up_residual_and_search_direction(const int &R_flag, const long int nodes, const long int &reserved_nodes, const Electric_para &electric_param, const vector<vector<double> > &KEFT, vector<double> &P, vector<double> &R, vector<double> &VEF)
{
    //Initialize the prescribed voltage boundary conditions
    //The magnitude of the voltage depends on the R_flag
    if (R_flag == 1) {
        
        //If R_flag is 1, then use real voltage (from the input parameters)
        if (!Get_voltage_vector(electric_param.applied_voltage, reserved_nodes, VEF)) {
            hout<<"Error in Set_up_residual_and_search_direction when calling Get_voltage_vector (R_falg="<<R_flag<<")"<< endl;
            return 0;
        }
    }
    else if (!R_flag) {
        
        //If R_flag is 0, then use the number of nodes as the magnitude for the voltage
        if (!Get_voltage_vector((double)nodes, reserved_nodes, VEF)) {
            hout<<"Error in Set_up_residual_and_search_direction when calling Get_voltage_vector (R_falg="<<R_flag<<")"<< endl;
            return 0;
        }
    }
    else {
        hout << "Error in Set_up_residual_and_search_direction. The R_flag has an invalid value: " << R_flag << endl;
        return 0;
    }
    hout << setwp(1,20) << "Maximum and minimum voltages = " << VEF.front() << ", " << VEF.back() << endl;
    
    //The residual vector is initialized with R = b - Ax0.
    //x0 is the initial guess. If we use x0 = 0 as initial guess then R = b
    //From the matrix equations b = - KEFT*VEF
    for (int i = 0; i < (int)KEFT.size(); i++) {
        
        for (int j = 1; j < (int)KEFT[i].size(); j++) {
            
            //The first element of VEF (i.e. VEF[0]) is zero,
            //so it can be skipped to save computational time
            R[i] = R[i] - (KEFT[i][j]*VEF[j]);
        }
        
        //The search direction of the CG is initialized with the initial value of the residual
        P[i] = R[i];
    }
    
    return 1;
}
//This function creates a voltage vector depending on the number of prescribed boundary conditios
int Direct_Electrifying::Get_voltage_vector(const double &volts, const long int &reserved_nodes, vector<double> &VEF)
{
    //Clear the vector of voltages
    for (int i = 0; i < reserved_nodes; i++) {
        VEF.push_back( volts*((double)i) );
    }
    
    return 1;
}
//This function solves the equation of the electric circuit using the
//Direct Electrifing Algorithm (DEA) and the Conjugate-Gradrient
int Direct_Electrifying::Solve_DEA_equations_CG_SSS(const long int &nodes, const long int &reserved_nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, vector<double> &P, vector<double> &R, vector<double> &VEF)
{
    //Preconditioner
    vector<double> M_inv(diagonal.size(),0.0);
    Jacobi_preconditioner(diagonal, M_inv);
    
    //Apply preconditoner
    vector<double> Y(R.size(),0.0);
    Apply_preconditioner(M_inv, R, P, Y);
    
    //Variables of the algorithm
    vector<double> AP; //Variable for the multiplication A*P
    AP.assign(nodes-reserved_nodes,0);
    double alpha, beta, rr0, rr;
    vector<double> voltages_sol(nodes-reserved_nodes, 0);
    
    //Maximum number of iterations for the CG
    long int max_iter = 10*nodes;
    //Iteration variable
    long int k;
    //Variable to check the status of the CG
    int test = 50000, test_inc = 50000;
    
    //Initial residual
    double R_dot_R = 0;
    if (!V_dot_v(R, R, R_dot_R)) {
        hout<<"Error in Solve_DEA_equations_CG_SSS when calling V_dot_v (0)"<<endl;
        return 0;
    }
    double R0 = 1.0E-10*sqrt(R_dot_R);
    //hout << "R0 = " << R0 << endl;
    
    //Preconditioned
    for (k = 1; k <= max_iter; k++) {
        
        //Calculate Ap
        spM_V_SSS(P, row_ptr, col_ind, diagonal, values, AP);
        
        //Calculate norm or residual of step k-1.
        //Will be used later as convergence criteria and to calculate beta
        if (!V_dot_v(Y, R, rr0)) {
            hout<<"Error in Solve_DEA_equations_CG_SSS when calling V_dot_v (1)"<<endl;
            return 0;
        }
        
        //Step length
        double dot_ = 0;
        if (!V_dot_v(Y, R, dot_)) {
            hout<<"Error in Solve_DEA_equations_CG_SSS when calling V_dot_v (2)"<<endl;
            return 0;
        }
        alpha = rr0/dot_;
        
        //Approximate solution
        //X = X + P*alpha;
        V_plus_aW(P, alpha, voltages_sol);
        
        //Residual
        //R = R - AP*alpha;
        V_plus_aW(AP, -alpha, R);
        
        //Calculate norm or residual of step k. Used as convergence criteria and to calculate beta
        if (!V_dot_v(R, R, rr)) {
            hout<<"Error in Solve_DEA_equations_CG_SSS when calling V_dot_v (3)"<<endl;
            return 0;
        }
        
        //Status update: print every test iterations
        if ( k == test){
            hout << "CG iteration " << k << endl;
            test = test + test_inc;
        }
        
        //Convergence criteria
        if (sqrt(rr) <= R0)
            break;
        
        //Update Y
        //Y = M_inv*R
        Componentwise_multiply(M_inv, R, Y);
        
        //Improvement of step
        if (!V_dot_v(Y, R, beta)) {
            hout<<"Error in Solve_DEA_equations_CG_SSS when calling V_dot_v (4)"<<endl;
            return 0;
        }
        beta = beta/rr0;
        
        //Search direction
        //P = Y + P*beta;
        W_plus_aV(Y, beta, P);
    }
    
    if (k >= max_iter)
        hout << "CG reached maximum number of iterations" << endl;
    hout << "CG iterations: " << k << endl;
    //hout << "RR = " << sqrt(rr) << endl;
    
    //Set up the voltages vector
    voltages.assign(nodes, 0);
    
    //Add the prescribed voltages
    for (long int i = 0; i < reserved_nodes; i++) {
        voltages[i] = VEF[i];
    }
    
    //Add the voltages from the solution
    for (long int i = reserved_nodes; i < nodes; i++) {
        voltages[i] = voltages_sol[i-reserved_nodes];
    }
    
    /*/Print the voltages
    Printer *Pr = new Printer;
    if (R_flag) {
        Pr->Print_1d_vec(voltages, "voltages_R.txt");
    } else {
        Pr->Print_1d_vec(voltages, "voltages_unit.txt");
    }
    delete Pr;//*/
    
    return 1;
}
//Calculate the Jacobi preconditioner M_inv
//M_inv[i] is simply 1/diagonal[i]
void Direct_Electrifying::Jacobi_preconditioner(const vector<double> &diagonal, vector<double> &M_inv)
{
    for (int i = 0; i < (int)diagonal.size(); i++) {
        //M_inv[i] is simply 1/diagonal[i]
        M_inv[i] = 1/diagonal[i];
    }
}
void Direct_Electrifying::Apply_preconditioner(const vector<double> &M_inv, const vector<double> &R, vector<double> &P, vector<double> &Y)
{
    //P = Y = M_inv*R
    for (int i = 0; i < (int)R.size(); i++) {
        P[i] = M_inv[i]*R[i];
        Y[i] = P[i];
    }
}
//This function performs a sparse-matrix-vector multiplication using the Symmetric Sparse Skyline (SSS) format
void Direct_Electrifying::spM_V_SSS(const vector<double> &V, const vector<long int> &rowptr, const vector<long int> &colind, const vector<double> &diagonal, const vector<double> &values, vector<double> &R)
{
    //Size of the system
    long int N = V.size();
    //Initialize result vector
    R.clear();
    R.assign(N,1);
    //SSS
    long int c;
    for (long int r = 0; r < N; r++) {
        R[r] = diagonal[r]*V[r];
        for (long int j = rowptr[r]; j < rowptr[r+1]; j++) {
            c = colind[j];
            R[r] = R[r] + values[j]*V[c];
            R[c] = R[c] + values[j]*V[r];
        }
    }
}
//This function calculates the dot product of two vectors
int Direct_Electrifying::V_dot_v(const vector<double> &A, const vector<double> &B, double &dot_product)
{
    if (A.size() != B.size()){
        hout << "Error in V_dot_v. Vectors A and B must have the same length. A.size()="<< A.size() << " B.size()=" << B.size() << endl;
        return 0;
    }
    
    //This variable will store the result of the dot product
    dot_product = 0;
    
    for (int i = 0; i < (int)A.size(); i++) {
        dot_product = dot_product + A[i]*B[i];
    }
    
    return 1;
}
//This function solves the following operation: V = V + a*W
//where V and W are vectors and a is a scalar
void Direct_Electrifying::V_plus_aW(const vector<double> &W, const double &a, vector<double> &V)
{
    for (int i = 0; i < (int)V.size(); i++) {
        V[i] = V[i] + a*W[i];
    }
}
//This function solves the following operation: V = W + a*V
//where V and W are vectors and a is a scalar
void Direct_Electrifying::W_plus_aV(const vector<double> &W, const double &a, vector<double> &V)
{
    for (int i = 0; i < (int)V.size(); i++) {
        V[i] = W[i] + a*V[i];
    }
}
//This function multiplies two vectors componentwise
void Direct_Electrifying::Componentwise_multiply(const vector<double> &vector_in1, const vector<double> &vector_in2, vector<double> &vector_out)
{
    for (int i = 0; i < (int)vector_in1.size(); i++) {
        vector_out[i] = vector_in1[i]*vector_in2[i];
    }
}
