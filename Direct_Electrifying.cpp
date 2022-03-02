//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Implementation of the Direct Electrifying Algorithm
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Direct_Electrifying.h"

//Calculate the voltage values at contact points and endpoints
int Direct_Electrifying::Compute_voltage_field(const int &n_cluster, const int &R_flag, const cuboid &window, const Simu_para &simu_para, const Electric_para &electric_param, const Cutoff_dist &cutoffs, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps)
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

    //hout << "global_nodes="<<global_nodes<<" reserved_nodes="<<reserved_nodes << endl;
    //Variables for using the SSS for the sparse matrix
    vector<long int> col_ind, row_ptr;
    vector<double> values, diagonal(global_nodes, 0);
    //Initialize the residual vector R
    vector<double> R(global_nodes-reserved_nodes,0);
    //Prescribed voltages applied to the sample
    vector<double> VEF(reserved_nodes, 0.0);
    
    //With the LM matrix, now fill the sparse stiffness matrix
    //Fill_sparse_stiffness_matrix (use R_flag)
    //hout<<"Fill_sparse_stiffness_matrix"<<endl;
    if (!Fill_sparse_stiffness_matrix(R_flag, global_nodes, reserved_nodes, cutoffs.van_der_Waals_dist, n_cluster, electric_param, HoKo, points_cnt, radii, structure_gnp, points_gnp, gnps, col_ind, row_ptr, values, diagonal, R, VEF)) {
        hout<<"Error in Compute_voltage_field when calling Fill_sparse_stiffness_matrix"<<endl;
        return 0;
    }
    
    /* /Check if actual resistances are used
    if (R_flag)
    {
        //Sort resistances
        sort(all_resistors.begin(), all_resistors.end());
        //Print file with sorted resistances
        Printer Pr; Pr.Print_1d_vec(all_resistors, "Resistors_" + to_string(n_cluster) + "_" + to_string(HoKo->family[n_cluster]) + ".txt");
    }// */
    
    //Create the approximation of the voltage vector when:
    //actual resistors are used
    //AND
    //the family is 0, 1, or 2
    vector<double> V_guess(global_nodes-reserved_nodes, 0.0);
    //Check if initial guess is needed
    if (R_flag && reserved_nodes == 2) {
        
        //Calculate an initial guess for the conjugate gradient
        //hout << "Initial_guess_for_CG" << endl;
        if (!Initial_guess_for_CG(n_cluster, reserved_nodes, HoKo->family[n_cluster], window, electric_param.applied_voltage, HoKo->clusters_cnt, HoKo->elements_cnt, points_cnt, HoKo->clusters_gnp, points_gnp, R, V_guess)) {
            hout<<"Error in Compute_voltage_field when calling Initial_guess_for_CG"<<endl;
            return 0;
        }

        //Update residual vector R
        //hout << "Update_residual_vector" << endl;
        if (!Update_residual_vector(col_ind, row_ptr, values, diagonal, V_guess, R)) {
            hout << "Error in Compute_voltage_field when calling Update_residual_vector" << endl;
            return 0;
        }
    }//
    
    //hout << "Solve_DEA_equations_CG_SSS"<<endl;
    //This is where the actual direct electrifying algorithm (DEA) takes place
    if (!Solve_DEA_equations_CG_SSS(global_nodes, reserved_nodes, simu_para.tolerance, col_ind, row_ptr, values, diagonal, R, VEF, V_guess)) {
        hout<<"Error in Compute_voltage_field when calling Solve_DEA_equations_CG_SSS"<<endl;
        return 0;
    }
    
    /* /Print the voltages
    Printer Pr;
    if (R_flag) {
        Pr.Print_1d_vec(voltages, "voltages_R_"+ to_string(HoKo->family[n_cluster])+ ".txt");
    } else {
        Pr.Print_1d_vec(voltages, "voltages_unit_" + to_string(HoKo->family[n_cluster]) + ".txt");
    }//*/
    
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
    //hout<<"Map_points_at_boundaries"<<endl;
    if (!Map_points_at_boundaries(HoKo->family[n_cluster], Cutwins->boundary_cnt_pts, LMM_cnts)) {
        hout<<"Error in LM_matrix_for_cnts when calling Map_points_at_boundaries (1)"<<endl;
        return 0;
    }
    
    //Iterate over all CNTs in the cluster
    for (int i = 0; i < (int)HoKo->clusters_cnt[n_cluster].size(); i++) {
        
        //Get the current CNT
        int CNT = HoKo->clusters_cnt[n_cluster][i];
        
        //Iterator for the elements in set
        set<long int>::iterator j = HoKo->elements_cnt[CNT].begin();
        //hout<<"HoKo->elements_cnt[CNT="<<CNT<<"].size="<<HoKo->elements_cnt[CNT].size()<<endl;
        //hout<<"P_begin="<<*j<<endl;
        
        //Iterator to find an element in LMM_cnts_boundary
        map<long int, long int>::iterator it;
        
        //Chek if the first point in the CNT element is at a boundary
        if (LMM_cnts.find(*j) == LMM_cnts.end()) {
            
            //Initial point of CNT is not at a boundary, so map it to a new node
            LMM_cnts[*j] = global_nodes;
            //hout<<"*j="<<*j<<" global_nodes="<<global_nodes<<endl;
            
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
        
        //hout<<"P_end="<<*j_end<<endl;
        //Chek if the last point in the CNT element is at a boundary
        if (LMM_cnts.find(*j_end) == LMM_cnts.end()) {
            
            //Last point of CNT is not at a boundary, so map it to a new node
            LMM_cnts[*j_end] = global_nodes;
            //hout<<"*j_end="<<*j_end<<" global_nodes="<<global_nodes<<endl;
            
            //Update the next available node
            global_nodes++;
        }
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
            //hout<<"P="<<P<<" node="<<n<<endl;
            
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
    //hout<<"Map_points_at_boundaries"<<endl;
    if (!Map_points_at_boundaries(HoKo->family[n_cluster], Cutwins->boundary_gnp_pts, LMM_gnps)) {
        hout<<"Error in LM_matrix_for_gnps when calling Map_points_at_boundaries"<<endl;
        return 0;
    }
    
    //Iterate over all GNPs in the cluster
    for (int i = 0; i < (int)HoKo->clusters_gnp[n_cluster].size(); i++) {
        
        //Get current GNP
        int gnp_i = HoKo->clusters_gnp[n_cluster][i];
        //hout << "structure_gnp[" << gnp_i << "].size()=" << structure_gnp[gnp_i].size() << endl;
        
        //Iterate over all points in GNP
        for (int j = 0; j < (int)structure_gnp[gnp_i].size(); j++) {
            
            //Get current point number
            long int P = structure_gnp[gnp_i][j];
            
            //Check it is not a boundary node with prescribed conditions
            if (LMM_gnps.find(P) == LMM_gnps.end()) {
                
                //P is not a boundary point
                
                //Add the mapping of the current GNP point to a node number
                LMM_gnps[P] = global_nodes;
                //hout << "LMM_gnps[P=" << P << "] = " << global_nodes << " gnp_i="<< gnp_i << endl;
                
                //Update the next available node
                global_nodes++;
            }
        }
    }
    
    return 1;
}
//This function creates the sparse stifness matrix that will be used to solve the sytem of equations
//The sparse version is more efficient computationally speaking
int Direct_Electrifying::Fill_sparse_stiffness_matrix(const int &R_flag, const long int &nodes, const long int &reserved_nodes, const double &d_vdw, const int &n_cluster, const Electric_para &electric_param, Hoshen_Kopelman *HoKo, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, vector<double> &R, vector<double> &VEF)
{
    //------------------------------------------------------------------------
    //Initially, generate a 2D vector that will be used to generate the 1D vectors for SSS
    //col_values[i] has the map col_j->value_ij, where i is the row number
    vector<map<long int, double>> col_values(nodes);
    
    //------------------------------------------------------------------------
    //Add contributions from particles
    
    //Fill the 2D matrices with the contributions of the CNTs when there are CNT clusters
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) 
    {
        //Add contributions from CNT resistors
        //hout << "Fill_2d_matrices_cnts"<<endl;
        if (!Fill_2d_matrices_cnts(R_flag, n_cluster, electric_param, HoKo, points_cnt, radii, col_values, diagonal)) 
        {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_cnts" << endl;
            return 0;
        }
        
        //Add contributions from CNT-CNT junctions if any
        if (HoKo->cluster_cnt_junctions.size() && HoKo->cluster_cnt_junctions[n_cluster].size()) 
        {
            //hout << "Fill_2d_matrices_cnt_junctions"<<endl;
            if (!Fill_2d_matrices_cnt_junctions(R_flag, reserved_nodes, d_vdw, electric_param, HoKo->cluster_cnt_junctions[n_cluster], HoKo->junctions_cnt, points_cnt, radii, LMM_cnts, col_values, diagonal)) 
            {
                hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_cnt_junctions" << endl;
                return 0;
            }
        }
    }
    
    //Set used to determine the nodes from CNT point in mixed junctions
    //map<long int, long int> points_cnt_rad;
    
    //Check if there are mixed junctions for n_cluster
    if (HoKo->cluster_mix_junctions.size() && HoKo->cluster_mix_junctions[n_cluster].size()) 
    {
        //Add contributions from mixed junctions
        //hout << "Fill_2d_matrices_mixed_junctions"<<endl;
        if (!Fill_2d_matrices_mixed_junctions(R_flag, reserved_nodes, d_vdw, electric_param, HoKo->cluster_mix_junctions[n_cluster], HoKo->junctions_mixed, points_cnt, radii, points_gnp, gnps, LMM_cnts, LMM_gnps, col_values, diagonal, points_cnt_rad)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_mixed_junctions" << endl;
            return 0;
        }
    }
    
    //Fill the 2D matrices with the contributions of the GNPs when there are GNP clusters
    if (HoKo->clusters_gnp.size() && HoKo->clusters_gnp[n_cluster].size()) 
    {
        //Add contributions from GNP resistors
        //hout<<"Fill_2d_matrices_gnp"<<endl;
        if (!Fill_2d_matrices_gnp(R_flag, electric_param, HoKo->clusters_gnp[n_cluster], points_gnp, structure_gnp, gnps, LMM_gnps, points_cnt_rad, col_values, diagonal)) 
        {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_gnp" << endl;
            return 0;
        }
        
        //Add contributions from GNP-GNP junctions, if any
        //hout<<"Fill_2d_matrices_gnp_junctions"<<endl;
        if (HoKo->cluster_gnp_junctions.size() && HoKo->cluster_gnp_junctions[n_cluster].size()) 
        {
            //hout << "junctions="<<HoKo->cluster_gnp_junctions[n_cluster].size() << endl;
            if (!Fill_2d_matrices_gnp_junctions(R_flag, d_vdw, electric_param, HoKo->cluster_gnp_junctions[n_cluster], HoKo->junctions_gnp, points_gnp, gnps, LMM_gnps, col_values, diagonal)) 
            {
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
    
    //hout<<"From_2d_to_1d_vectors"<<endl;
    //Set up variables for the Conjugate Gradient Algorithm:
    //1D vectors for SSS format: col_ind, row_ptr, values
    //Search direction (P) and residual vector (R)
    if (!From_2d_to_1d_vectors(reserved_nodes, nodes, R_flag, electric_param, col_values, col_ind, row_ptr, values, diagonal, R, VEF)) {
        hout<<"Error in Fill_sparse_stiffness_matrix when calling From_2d_to_1d_vectors"<<endl;
        return 0;
    }
    
    /* /Print some vectors
    Printer Pr;
    string str = "-" + to_string(R_flag) + "-" + to_string(HoKo->family[n_cluster]);
    Pr.Print_1d_vec(diagonal, "diagonal" + str + ".txt");
    Pr.Print_1d_vec(row_ptr, "row_ptr" + str + ".txt");
    Pr.Print_1d_vec(col_ind, "col_ind" + str + ".txt");
    Pr.Print_1d_vec(values, "values" + str + ".txt");
    Pr.Print_1d_vec(R, "R" + str + ".txt");// */

    /* /Print some vectors with sorted values of the non-zero entries of the stifness matrix
    Printer Pr;
    string str = "-" + to_string(R_flag) + "-" + to_string(HoKo->family[n_cluster]);
    vector<double> diag_sorted = diagonal;
    sort(diag_sorted.begin(), diag_sorted.end());
    Pr.Print_1d_vec(diag_sorted, "diag_sorted" + str + ".txt");
    vector<double> val_sorted = values;
    sort(val_sorted.begin(), val_sorted.end());
    Pr.Print_1d_vec(val_sorted, "val_sorted" + str + ".txt");// */

    //Check there are no zeros in the diagonal
    //At this point, the diagonal vector is in its final form
    vector<int> zero_idx;
    vector<double> zero_vals;
    for (size_t i = 0; i < diagonal.size(); i++) {
        if (abs(diagonal[i]) < Zero) {
            //Add the reserved nodes, since at this points those have been removed
            //from the diagonal vector
            zero_idx.push_back(i + reserved_nodes);
            zero_vals.push_back(diagonal[i]);
        }
    }
    
    if (zero_idx.size()) {
        hout << "Error in Fill_sparse_stiffness_matrix. There are zero elements in the diagonal of the sitffness matrix. \nThis means that there are nodes not connected no any other node." << endl;
        hout << "Nodes with zero on the diagonal: ";
        for (size_t i = 0; i < zero_idx.size(); i++) {
            hout << zero_idx[i] << " ";
        }
        hout << endl;
        hout << "Zeros on the diagonal: ";
        for (size_t i = 0; i < zero_vals.size(); i++) {
            hout << zero_vals[i] << " ";
        }
        hout << endl;
        return 0;
    }

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
        //hout<<"P1="<<P1<<" node="<<LMM_cnts[P1]<<endl;
        //hout<<"CNT="<<CNT<<endl;
        
        //Check for the case of three points in the element
        if (HoKo->elements_cnt[CNT].size() == 3 &&
            LMM_cnts[P1] == LMM_cnts[*(HoKo->elements_cnt[CNT].rbegin())]) {
            
            //Go to the case of three points with the endpoints in the same boundary
            //hout<<"Three_points_in_element"<<endl;
            if (!Three_points_in_element(R_flag, electric_param, HoKo->elements_cnt[CNT], points_cnt, radii[CNT], col_values, diagonal)) {
                hout<<"Error in Fill_2d_matrices_cnts when calling Three_points_in_element"<<endl;
                return 0;
            }
        }
        else {
            
            //Iterate over the points in the elements vector of CNT
            for (it++; it != HoKo->elements_cnt[CNT].end(); it++) {
                
                //Get the second point of the element
                long int P2 = *it;
                
                //Calculate resistance depending of the R_flag
                double Re;
                if (!Calculate_resistance_cnt(R_flag, points_cnt, P1, P2, radii[CNT], electric_param, Re)) {
                    hout<<"Error in Fill_2d_matrices_cnts when calling Calculate_resistance_cnt"<<endl;
                    return 0;
                }
                
                //Get the nodes of the points
                long int node1 = LMM_cnts[P1];
                long int node2 = LMM_cnts[P2];
                //hout<<"Pi="<<P2<<" node="<<node2<<endl;
                
                //Calculate inverse of resistance if needed
                double Re_inv = (R_flag)? 1/Re : Re;
                
                //Add the resistance value to the corresponding 2D vectors
                if(!Add_new_elements_to_2d_sparse_matrix(node1, node2, Re_inv, col_values, diagonal)) {
                    hout<<"Error in Fill_2d_matrices_cnts when calling Add_elements_to_2d_sparse_matrix"<<endl;
                    return 0;
                }
                
                //Set the second point as the first point for the next iteration
                P1 = P2;
            }
        }
    }
    
    return 1;
}
//This function adds resistors to the stiffness matrix for the case when there are three
//points in the CNT element and the two endpoints are connected to the same boundary
int Direct_Electrifying::Three_points_in_element(const int &R_flag, const Electric_para &electric_param, const set<long int> &cnt_element, const vector<Point_3D> &points_cnt, const double &radius, vector<map<long int, double>> &col_values, vector<double> &diagonal)
{
    //Get the three points
    set<long int>::const_iterator it = cnt_element.begin();
    long int P1 = *it; ++it;
    long int P2 = *it; ++it;
    long int P3 = *it;
    
    //Calculate resistance of segment 1
    double Re1_inv = 1.0;
    if (R_flag == 1) {
        
        //hout<<"Calculate_resistance_cnt 1"<<endl;
        if (!Calculate_resistance_cnt(R_flag, points_cnt, P1, P2, radius, electric_param, Re1_inv)) {
            hout<<"Error in Three_points_in_element when calling Calculate_resistance_cnt Re1"<<endl;
            return 0;
        }
        
        //Calculate the inverse of the resistance, as this value is the one needed
        //for the stiffness matrix
        Re1_inv = 1/Re1_inv;
    }
    
    //Calculate resistance of segment 2
    double Re2_inv = 1.0;
    if (R_flag == 1) {
        
        //hout<<"Calculate_resistance_cnt 2"<<endl;
        if (!Calculate_resistance_cnt(R_flag, points_cnt, P2, P3, radius, electric_param, Re2_inv)) {
            hout<<"Error in Three_points_in_element when calling Calculate_resistance_cnt Re2"<<endl;
            return 0;
        }
        
        //Calculate the inverse of the resistance, as this value is the one needed
        //for the stiffness matrix
        Re2_inv = 1/Re2_inv;
    }
    
    //Get the nodes of the first two points
    //The third node is not needed as it is the same as the first node
    long int node1 = LMM_cnts[P1];
    long int node2 = LMM_cnts[P2];
    //hout<<"P1="<<P1<<" P2="<<P2<<" P3="<<P3<<" CNT="<<points_cnt[P2].flag<<endl;
    //hout<<"node1="<<node1<<" node2="<<node2<<" node3="<<LMM_cnts[P3]<<endl;
    
    //Add contributions to the diagonal of the stiffness matrix
    diagonal[node1] += Re1_inv + Re2_inv;
    diagonal[node2] += Re1_inv + Re2_inv;
    
    //Add contributions to the off diagonal entries of the stiffness matrix
    col_values[node2][node1] = -(Re1_inv + Re2_inv);
    
    return 1;
}
//This function calculates the length of a CNT segment that corresponds to one element (resistor)
//Using the resisitivity as input parameter the resistance of the CNT segment is calculated
int Direct_Electrifying::Calculate_resistance_cnt(const int &R_flag, const vector<Point_3D> &points_cnt, const long int &P1, const long int &P2, const double &radius, const Electric_para& electric_param, double &Re)
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
        Re = electric_param.resistivity_CNT*length/(PI*radius*radius);

        //Scaling factor
        Re = Re * electric_param.scaling_R;

        //Add resistance to vector of all resistors
        //all_resistors.push_back(Re);
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
//This function adds new elements to the 2D sparse stiffness matrix
int Direct_Electrifying::Add_new_elements_to_2d_sparse_matrix(const long int &node1, const long int &node2, const double &Re_inv, vector<map<long int, double> > &col_values, vector<double> &diagonal)
{
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
//This function adds elements to the 2D sparse stiffness matrix when an element connecting the
//same nodes has already been added
int Direct_Electrifying::Add_to_existing_elements_in_2d_sparse_matrix(const long int &node1, const long int &node2, const double &Re_inv, vector<map<long int, double> > &col_values, vector<double> &diagonal)
{
    //Add the diagonal elements of the stiffness matrix
    diagonal[node1] += Re_inv;
    diagonal[node2] += Re_inv;
    
    //Add the off diagonal elements of the stiffness matrix
    if (node1 > node2) {
        
        //node2 is the column index and node1 is the row index
        //In row node1, add to the map node2-> -Re_inv
        col_values[node1][node2] = col_values[node1][node2] -Re_inv;
        
    } else {
        
        //node1 is the column index and node2 is the row index
        //In row node2, add to the map node1-> -Re_inv
        col_values[node2][node1] = col_values[node2][node1] -Re_inv;
    }
    
    return 1;
}
//This function checks if an element was already added to the 2D sparse stiffness matrix
//If not already in the matrix, it is added
//If already in the matrix, the contribution of the new element is added to the existing value
int Direct_Electrifying::Add_new_or_to_existing_elements_in_2d_sparse_matrix(const long int &node1, const long int &node2, const double &Re_inv, vector<map<long int, double> > &col_values, vector<double> &diagonal)
{
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
int Direct_Electrifying::Fill_2d_matrices_cnt_junctions(const int &R_flag, const long int &reserved_nodes, const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_cnt_junctions_i, const vector<Junction> &junctions_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const map<long int, long int> &LMM_cnts, vector<map<long int, double> > &col_values, vector<double> &diagonal)
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
        //hout<<"P1="<<P1<<" CNT1="<<CNT1<<" P2="<<P2<<" CNT2="<<CNT2<<endl;
        //hout<<"idx="<<idx<<" junctions_cnt.size="<<junctions_cnt.size()<<endl;
        
        //Calculate the junction resistance
        double Re_inv = 1.0;
        if (R_flag == 1) {
            if (!Calculate_junction_resistance(junctions_cnt[idx], radii[CNT1], radii[CNT2], electric_param, Re_inv)) {
                hout<<"Error in Fill_2d_matrices_cnt_junctions when calling Calculate_junction_resistance"<<endl;
                return 0;
            }

            //Add resistance to vector of all resistors
            //all_resistors.push_back(Re_inv);
            
            //Calculate inverse of resistance
            Re_inv = 1/Re_inv;
        }
        else if (R_flag != 0) {
            hout << "Error in Fill_2d_matrices_cnt_junctions. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
            return 0;
        }
        
        //Get node numbers
        long int node1 = LMM_cnts.at(P1);
        long int node2 = LMM_cnts.at(P2);
        //hout<<"node1="<<node1<<" node2="<<node2<<endl;
        /* /Check for small numbers
        if (Re_inv < 1e-13) {
            hout << "R_inv=" << Re_inv << " d=" << junctions_cnt[idx].junction_dist << " d_P=" << points_cnt[P1].distance_to(points_cnt[P2]) << endl;
            hout << "P1=" << P1 << " CNT1=" << CNT1 << " P2=" << P2 << " CNT2=" << CNT2 << endl;
            hout << "r_cnt1=" << radii[CNT1] << " r_cnt2=" << radii[CNT2] << endl;
            hout << "node1=" << node1 << " node2=" << node2 << endl;
            double Re = electric_param.C1 * junctions_cnt[idx].junction_dist * exp(electric_param.C2 * junctions_cnt[idx].junction_dist) / (radii[CNT1]*radii[CNT2]*4.0);
            hout << "Re=" << Re << " 1/Re=" << 1.0 / Re << endl;
        }// */
        
        //Check if any of the nodes is at a boundary with prescribed voltage
        if (node1 >= reserved_nodes && node2 >= reserved_nodes) {
            
            //Add junction resistance to sparse stiffness matrix
            //hout<<"Add_new_elements_to_2d_sparse_matrix i="<<i<<endl;
            if (!Add_new_elements_to_2d_sparse_matrix(node1, node2, Re_inv, col_values, diagonal)) {
                hout<<"Error in Fill_2d_matrices_cnt_junctions when calling Add_elements_to_2d_sparse_matrix"<<endl;
                return 0;
            }
        }
        //Add the junction resistors when calculating the voltage field and there are more than 2
        //boundaries with prescribed voltage
        //Otherwise, the junction is ignored
        else if (reserved_nodes > 2) {
            
            //Add junction resistance to existing elements in sparse stiffness matrix
            //hout<<"Add_new_elements_to_2d_sparse_matrix"<<endl;
            if (!Add_to_existing_elements_in_2d_sparse_matrix(node1, node2, Re_inv, col_values, diagonal)) {
                hout<<"Error in Fill_2d_matrices_cnt_junctions when calling Add_to_existing_elements_in_2d_sparse_matrix"<<endl;
                return 0;
            }
        }
    }
    
    return 1;
}
//This function calculates the junction resistance between either two CNTs, two GNPs, or
//a CNT and a GNP
int Direct_Electrifying::Calculate_junction_resistance(const Junction &j, const double &l1, const double &l2, const struct Electric_para &electric_param, double &Re)
{
    //Check if a constant resistance is used
    if (electric_param.junction_type == "constant") 
    {
        //Take the constant value from the electrical parameters
        Re = electric_param.junction_resistance;
    }
    else if (electric_param.junction_type == "exponential") 
    {

        //Variables to store the CNT diameter or GNP thickness, which is neede to calculate
        //the area of the junction
        double jl1 = l1, jl2 = l2;
        
        //Check if the first particle of the junction is a CNT
        if (j.type1 == "CNT") 
        {
            //Get the CNT diameter
            jl1 = jl1 + l1;
        }
        
        //Check if the second particle of the junction is a CNT
        if (j.type2 == "CNT")
        {
            //Get the CNT diameter
            jl2 = jl2 + l2;
        }

        //Junction area is considered to be a square
        //At this point in jl1 and jl2 the CNT diameter of GNP thickness is stored
        //So calculate the area of the junction assuming a square using the smallest
        //CNT diameter or GN thickness, which will the one limiting the flow of electrons
        double A = (jl1 - jl2 < Zero)? jl1 * jl1 : jl2 * jl2;
        
        //Calculate junction resistance using the pre-calculated constants
        Re = electric_param.C1 * j.junction_dist * exp(electric_param.C2*j.junction_dist)/A;
    }
    else {
        hout<<"Error in Calculate_junction_resistance: invalid junction type. Input was: "<<electric_param.junction_type<<endl;
        return 0;
    }

    //Scaling factor
    Re = Re * electric_param.scaling_R;
    
    return 1;
}
//This function adds the contributions of the junctions between a CNT and a GNP
int Direct_Electrifying::Fill_2d_matrices_mixed_junctions(const int &R_flag, const long int &reserved_nodes, const double &d_vdw, const Electric_para &electric_param, const vector<int> cluster_mix_junctions_i, const vector<Junction> &junctions_mixed, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const map<long int, long int> &LMM_cnts, const map<long int, long int> &LMM_gnps, vector<map<long int, double> > &col_values, vector<double> &diagonal, map<long int, double> &points_cnt_rad)
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
        double Re_inv = 1.0;
        if (R_flag == 1) {
            if (!Calculate_junction_resistance(junctions_mixed[idx], radii[cnt_n], gnps[gnp_n].t, electric_param, Re_inv)) {
                hout<<"Error in Fill_2d_matrices_mixed_junctions when calling Calculate_junction_resistance"<<endl;
                return 0;
            }

            //Add resistance to vector of all resistors
            //all_resistors.push_back(Re_inv);
            
            //Calculate inverse of resistance
            Re_inv = 1/Re_inv;
        }
        else if (R_flag != 0) {
            hout << "Error in Fill_2d_matrices_mixed_junctions. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
            return 0;
        }
        
        //Get node numbers
        long int node1 = LMM_cnts.at(Pcnt);
        long int node2 = LMM_gnps.at(Pgnp);
        /* /Check for small numbers
        if (Re_inv < 1e-13) {
            hout << "R_inv=" << Re_inv << " d=" << junctions_mixed[idx].junction_dist << " d_P=" << points_gnp[Pgnp].distance_to(points_cnt[Pcnt]) << endl;
            hout << "Pcnt=" << Pcnt << " cnt_n=" << cnt_n << " Pgnp=" << Pgnp << " gnp_n=" << gnp_n << endl;
            hout << "r_cnt=" << radii[cnt_n] << endl;
            hout << "l_gnp=" << gnps[gnp_n].l << " t_gnp=" << gnps[gnp_n].t << endl;
            hout << "node1=" << node1 << " node2=" << node2 << endl;
            double Re = electric_param.C1 * junctions_mixed[idx].junction_dist * exp(electric_param.C2 * junctions_mixed[idx].junction_dist) / (gnps[gnp_n].t * radii[cnt_n] * 2.0);
            hout << "Re=" << Re << " 1/Re=" << 1.0 / Re << endl;
        }// */
        
        //Check if any of the nodes is at a boundary with prescribed voltage
        if (node1 >= reserved_nodes && node2 >= reserved_nodes) {
            
            //Add junction resistance to sparse stiffness matrix
            //hout<<"Add_new_elements_to_2d_sparse_matrix i="<<i<<endl;
            if (!Add_new_elements_to_2d_sparse_matrix(node1, node2, Re_inv, col_values, diagonal)) {
                hout<<"Error in Fill_2d_matrices_mixed_junctions when calling Add_elements_to_2d_sparse_matrix"<<endl;
                return 0;
            }
        }
        //Add the junction resistors when calculating the voltage field and there are more than 2
        //boundaries with prescribed voltage
        //Otherwise, the junction is ignored
        else if (reserved_nodes > 2) {
            
            //Add junction resistance to existing elements in sparse stiffness matrix
            //hout<<"Add_new_elements_to_2d_sparse_matrix"<<endl;
            if (!Add_to_existing_elements_in_2d_sparse_matrix(node1, node2, Re_inv, col_values, diagonal)) {
                hout<<"Error in Fill_2d_matrices_mixed_junctions when calling Add_to_existing_elements_in_2d_sparse_matrix"<<endl;
                return 0;
            }
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
        
        //current GNP
        int gnp_i = cluster_gnp[i];
        //hout<<endl<<"gnp_i="<<gnp_i<<" structure_gnp[gnp_i].size="<< structure_gnp[gnp_i].size() << " gnps[gnp_i].triangulation.size()="<< gnps[gnp_i].triangulation.size()<<endl;
        //for (size_t i = 0; i < structure_gnp[gnp_i].size(); i++) hout << "structure_gnp[gnp_i="<< gnp_i<<"][" << i << "]=" << structure_gnp[gnp_i][i] << endl;
        
        //Check if the triangulation needs to be performed
        //Triangulation is performed when calculating the backbone (R_flag == 0)
        //or when calculating the electrical resistance and there is already a triangulation
        if (!R_flag || gnps[gnp_i].triangulation.empty()) {
            
            //Perform triangulation
            if (!dt.Generate_3d_trangulation(points_gnp, structure_gnp[gnp_i], gnps[gnp_i])) {
                hout << "Error in Fill_2d_matrices_gnp when calling delaunay->Generate_3d_trangulation" << endl;
                return 0;
            }
            //hout << "Triangulation " << i << ", " << structure_gnp[gnp_i].size() << " points in GNP " << gnp_i<<", " << gnps[gnp_i].triangulation.size() << " edges"<< endl;
            if (structure_gnp[gnp_i].size() > 4 && !gnps[gnp_i].triangulation.size())
            {
                hout << "Error in Fill_2d_matrices_gnp when calling Generate_3d_trangulation. Bowyer-Watson algorithm is not generating a triangulation for more than 4 points." << endl;
                return 0;
            }
        }
        
        //Add elements from the triangulation
        //Scan triangulation edges in reverse order since some edges might need to be deleted
        for (int j = (int)gnps[gnp_i].triangulation.size() - 1; j >= 0 ; j--) {

            //hout << "Edge " << j << endl;
            
            //Get the two vertices of the triangulation
            long int v1 = gnps[gnp_i].triangulation[j].v1;
            long int v2 = gnps[gnp_i].triangulation[j].v2;
            
            //Get the node numbers
            long int node1 = LMM_gnps.at(v1);
            //hout<<"v1="<<v1<<" node1="<<node1<<endl;
            long int node2 = LMM_gnps.at(v2);
            //hout<<"v2="<<v2<<" node2="<<node2<<endl;
            
            //Initialize resistor with unit resistance
            double Re_inv = 1.0;
            
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
                
                //hout<<"rad1="<<rad1<<" rad2="<<rad2<<endl;
                //hout<<"points_gnp[v1]="<<points_gnp[v1].str()<<" points_gnp[v2]="<<points_gnp[v2].str()<<endl;
                //Calculate the triangulation resistor, i.e., the resistance of
                //the "conduction band" in the GNP
                if (!Calculate_resistance_gnp(points_gnp[v1], points_gnp[v2], rad1, rad2, electric_param, Re_inv)) {
                    hout << "GNP_i=" << gnp_i << endl;
                    hout << "Error in Fill_2d_matrices_gnp when calling Calculate_resistance_gnp" << endl;
                    return 0;
                }
                
                //Add resistance to vector of all resistors
                //all_resistors.push_back(Re_inv);

                //Calculate inverse of resistance
                //hout << "R=" << Re_inv << " ";
                Re_inv = 1/Re_inv;
                //hout << "R_inv=" << Re_inv << endl;
            }
            else if (R_flag != 0) {
                hout << "Error in Fill_2d_matrices_gnp. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
                return 0;
            }
            
            //hout << "Add_new_elements_to_2d_sparse_matrix" << endl;
            if (!Add_new_elements_to_2d_sparse_matrix(node1, node2, Re_inv, col_values, diagonal)) {
                hout << "Error in Fill_2d_matrices_gnp when calling Add_elements_to_2d_sparse_matrix" << endl;
                return 0;
            }
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

    if (L < Zero)
    {
        hout << "Error in Calculate_resistance_gnp. GNP points to calculate GNP resistor are too close or the same." << endl;
        hout<<"This results in a resistor with zero length." << endl;
        hout << "P1=" << P1.str() << endl;
        hout << "P2=" << P2.str() << endl;
        hout << "P1.distance_to(P2)=" << L << endl;
        return 0;
    }
    
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

    //Scaling factor
    Re = Re * electric_param.scaling_R;
    
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
        //hout << "Junction: GNP1=" << GNP1 << " GNP2=" << GNP2 << endl;
        
        //Calculate the junction resistance
        double Re_inv = 1.0;
        if (R_flag == 1) {
            if (!Calculate_junction_resistance(junctions_gnp[idx], gnps[GNP1].t, gnps[GNP2].t, electric_param, Re_inv)) {
                hout<<"Error in Fill_2d_matrices_gnp_junctions when calling Calculate_junction_resistance"<<endl;
                return 0;
            }

            //Add resistance to vector of all resistors
            //all_resistors.push_back(Re_inv);
            
            //Calculate inverse of resistance
            //hout << "Rj=" << Re_inv << " ";
            Re_inv = 1/Re_inv;
            //hout << "Rj_inv=" << Re_inv << endl;
        }
        else if (R_flag != 0) {
            hout << "Error in Fill_2d_matrices_gnp_junctions. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
            return 0;
        }
        
        //Get node numbers
        long int node1 = LMM_gnps.at(P1);
        long int node2 = LMM_gnps.at(P2);
        /* /Check for small numbers
        if (Re_inv < 1e-13) {
            hout << "R_inv=" << Re_inv << " d=" << junctions_gnp[idx].junction_dist <<" d_P="<< points_gnp[P1].distance_to(points_gnp[P2]) << endl;
            hout << "GNP1=" << GNP1 << " GNP2=" << GNP2 << endl;
            hout << "l_GNP1=" << gnps[GNP1].l << " t_GNP1=" << gnps[GNP1].t << endl;
            hout << "l_GNP2=" << gnps[GNP2].l << " t_GNP2=" << gnps[GNP2].t << endl;
            hout << "node1=" << node1 << " node2=" << node2 << endl;
            double Re = electric_param.C1 * junctions_gnp[idx].junction_dist * exp(electric_param.C2 * junctions_gnp[idx].junction_dist) / (gnps[GNP1].t* gnps[GNP1].t);
            hout << "Re=" << Re << " 1/Re=" << 1.0 / Re << endl;
        }// */
        
        //Add junction resistance to sparse stiffness matrix
        if (!Add_new_elements_to_2d_sparse_matrix(node1, node2, Re_inv, col_values, diagonal)) {
            hout<<"Error in Fill_2d_matrices_mixed_junctions when calling Add_elements_to_2d_sparse_matrix"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function transforms the 2D vector (vector of maps) col_values into 1D vectors
//so they can be in the SSS format and make the matrix-vector multiplications faster
//For reference, the 2D stiffness matrix has this form:
//
//   | KE    KEF |
//   | KEFT  KF  |
//
//For the CG algorithm, this function sets up the search direction (P) and
//residual vector (R), and also the vector of prescribed voltages VEF
int Direct_Electrifying::From_2d_to_1d_vectors(const long int &reserved_nodes, const long int &nodes, const int &R_flag, const Electric_para &electric_param, const vector<map<long int, double> > &col_values, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, vector<double> &R, vector<double> &VEF)
{
    //Initialize the prescribed voltage boundary conditions
    //hout<<"Get_voltage_vector"<<endl;
    if (!Get_voltage_vector(reserved_nodes, nodes, R_flag, electric_param, VEF)) {
        hout<<"Error in From_2d_to_1d_vectors when calling Get_voltage_vector"<<endl;
        return 0;
    }
    
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
                
                //When the column index is less than reserved_nodes, then the value of Kij
                //actually corresponds to the matrix KEFT
                //The column index is (of course) the column index of KEFT as
                //in this else-statement we have that 0<=col<reserved_nodes
                //The row is the index iterator i-reserved_nodes
                //Thus we have that for the KEFT matrix:
                //KEFT[i-reserved_nodes][col] = Kij;
                
                //The residual vector is initialized with R = b - Ax0,
                //where x0 is the initial guess. If we use x0 = 0 as initial guess then R = b
                //From the matrix equations b = - KEFT*VEF
                //Thus, the KEFT matrix is used to initialize the residual vector as follows:
                //hout << "\tR[" << i - reserved_nodes << "]=" << R[i - reserved_nodes] << " - " << Kij << "*" << VEF[col] << endl;
                R[i-reserved_nodes] = R[i-reserved_nodes] - (Kij*VEF[col]);
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
//This function calculates the initial guess for the Conjugate Gradien (CG)
int Direct_Electrifying::Initial_guess_for_CG(const int &n_cluster, const long int &reserved_nodes, const int &family, const cuboid &window_geom, const double &V_app, const vector<vector<int> > &clusters_cnt, const vector<set<long int> > &elements_cnt, const vector<Point_3D> &points_cnt, const vector<vector<int> > &clusters_gnp, const vector<Point_3D> &points_gnp, vector<double> &R, vector<double> &V_guess)
{
    
    //Select the length of the window along the direction indicated by the family
    double L_window = Select_coordinate(family, Point_3D(window_geom.len_x, window_geom.wid_y, window_geom.hei_z));
    
    //Select the maximum coordinate of the window along the direction indicated by the family
    double max_coord = Select_coordinate(family, Point_3D(window_geom.max_x, window_geom.max_y, window_geom.max_z));
    
    //Precalculate V_app/L_window
    double V_P_coord = V_app/L_window;

    //Add the initial guess for the nodes tha correspond to CNT points (if any)
    if (clusters_cnt.size() && clusters_cnt[n_cluster].size()) {
        hout << "Initial_guess_for_cnt_nodes" << endl;
        if (!Initial_guess_for_cnt_nodes(n_cluster, reserved_nodes, family, V_P_coord, max_coord, clusters_cnt, elements_cnt, points_cnt, V_guess))
        {
            hout << "Error in Initial_guess_for_CG when calling Initial_guess_for_cnt_nodes" << endl;
            return 0;
        }
    }
    
    //Print the initial guess for the voltages
    /*Printer Pr;
    string filename = "V_guess_" + to_string(family) + ".txt";
    Pr.Print_1d_vec(V_guess, filename);*/

    return 1;
}
//This function selects the coordinate corresponding to a family:
//family 0 = x-coordinate
//family 1 = y-coordinate
//family 2 = z-coordinate
double Direct_Electrifying::Select_coordinate(const double& family, const Point_3D& P)
{
    return ( (family==0)? P.x: (family==1)? P.y: P.z);
}
//Calculate the initial guess for the nodes that correspond to CNT points
int Direct_Electrifying::Initial_guess_for_cnt_nodes(const int& n_cluster, const long int& reserved_nodes, const int& family, const double& V_P_coord, const double& max_coord, const vector<vector<int> >& clusters_cnt, const vector<set<long int> >& elements_cnt, const vector<Point_3D>& points_cnt, vector<double>& V_guess)
{
    //Scan every CNT in the cluster
    for (int i = 0; i < (long int)clusters_cnt[n_cluster].size(); i++) {

        //Current CNT in cluster
        int CNT = clusters_cnt[n_cluster][i];

        //Iterate over the points in the elements vector of CNT
        for (set<long int>::const_iterator it = elements_cnt[CNT].begin();
            it != elements_cnt[CNT].end(); it++) {

            //Get the current point of the element
            long int Pi = *it;

            //Get the node number of the point
            long int node_i = LMM_cnts.at(Pi);

            //Check if node_i is not a reserved node
            if (node_i >= 2) {

                //Select the proper coordinate
                double P_coord = Select_coordinate(family, points_cnt[Pi]);

                //Set the voltage at node_i porportional to the coordinate P_coord
                //The vector used in CG does not have the reserved nodes, thus these
                //are subtracted from the index
                V_guess[node_i - reserved_nodes] = V_P_coord * (max_coord - P_coord);

            }
            //Voltages for nodes 0 and 1 have already been assigned, so there is nothing
            //to do if those nodes are found in this function
        }
    }

    return 1;
}
//This function creates a voltage vector depending on the number of prescribed boundary conditios
int Direct_Electrifying::Get_voltage_vector(const long int &reserved_nodes, const long int &nodes, const int &R_flag, const Electric_para &electric_param, vector<double> &VEF)
{
    //Variable to store the magnitude of the voltage
    double V;
    
    //The magnitude of the voltage depends on the R_flag
    if (R_flag == 1) {
        
        //If R_flag is 1, then use real voltage (from the input parameters)
        V = electric_param.applied_voltage;
    }
    else if (!R_flag) {
        
        //If R_flag is 0, then use the number of nodes as the magnitude for the voltage
        V = (double)nodes;
    }
    else {
        hout << "Error in Get_voltage_vector. The R_flag has an invalid value: " << R_flag << endl;
        return 0;
    }
    
    //Fill the vector of prescribed voltages, ignore the first entry as it is already 0
    for (int i = 1; i < reserved_nodes; i++) {
        VEF[i] = VEF[i-1] + V;
    }
    hout << setwp(1,20) << "Maximum and minimum voltages = " << VEF.front() << ", " << VEF.back() << endl;
    
    return 1;
}
//After calculating an initial guess, this function updates the residual vector R
int Direct_Electrifying::Update_residual_vector(const vector<long int>& col_ind, const vector<long int>& row_ptr, const vector<double>& values, const vector<double>& diagonal, const vector<double>& V_guess, vector<double>& R)
{
    //Calculate A*x0 needed to calculate the residual
    //The matrix-vector multiplication is done using the SSS format

    //Size of the system
    long int N = (long int)V_guess.size();

    //Initialize result vector
    vector<double> Ax0(N, 1.0);
    //SSS
    //hout<<"SSS, N="<<N<<" col_ind.size="<<col_ind.size()<<" row_ptr.size="<<row_ptr.size()<<endl;
    for (long int r = 0; r < N; r++) {

        Ax0[r] = diagonal[r] * V_guess[r];
        //hout<<"r="<<r<<endl;
        for (long int j = row_ptr[r]; j < row_ptr[r + 1]; j++) {

            //Get the column index c
            long int c = col_ind[j];
            //hout<<" c="<<c<<endl;
            Ax0[r] = Ax0[r] + values[j] * V_guess[c];
            Ax0[c] = Ax0[c] + values[j] * V_guess[r];
        }
    }

    //Calculate the residual R = b - A*x0
    //Variable R already has stored the value of b so I just need to do: R <- b - A*x0
    for (long int i = 0; i < N; i++) {
        R[i] = R[i] - Ax0[i];
    }

    return 1;
}
//This function solves the equation of the electric circuit using the
//Direct Electrifing Algorithm (DEA) and the Conjugate-Gradrient
int Direct_Electrifying::Solve_DEA_equations_CG_SSS(const long int &nodes, const long int &reserved_nodes, const double &tolerance, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, vector<double> &R, vector<double> &VEF, vector<double> &voltages_sol)
{
    //Preconditioner
    vector<double> M_inv(diagonal.size(),0.0);
    Jacobi_preconditioner(diagonal, M_inv);
    //Printer Pr;
    //Pr.Print_1d_vec(M_inv, "M_inv.txt");
    
    //Apply preconditoner
    vector<double> Y(R.size(),0.0), P(R.size(),0.0);
    Apply_preconditioner(M_inv, R, P, Y);
    //Pr.Print_1d_vec(P, "P_prec.txt");
    //Pr.Print_1d_vec(Y, "Y_prec.txt");
    
    //Variables of the algorithm
    vector<double> AP; //Variable for the multiplication A*P
    AP.assign(nodes-reserved_nodes,0);
    double alpha, beta, rr0, rr;
    //vector<double> voltages_sol(nodes-reserved_nodes, 0);
    
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
    double R0 = tolerance*sqrt(R_dot_R);
    //hout<<"R0="<<R0<<" tolerance="<<tolerance<< endl;
    
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
        //hout<<"k= "<<k<<"\trr0="<<rr0<<endl;
        
        //Step length
        double dot_ = 0;
        if (!V_dot_v(P, AP, dot_)) {
            hout<<"Error in Solve_DEA_equations_CG_SSS when calling V_dot_v (2)"<<endl;
            return 0;
        }
        //hout<<"\tdot_="<<dot_<<endl;
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
        //hout << "\trr=" << rr << endl;
        
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
        //hout << "\tbeta=" << beta;
        beta = beta/rr0;
        //hout << " beta/rr0=" << beta << endl;
        
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
