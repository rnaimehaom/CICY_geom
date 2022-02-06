//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Determine the backbone network and dead branches in the percolation network
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Backbone_Network.h"


int Backbone_Network::Determine_backbone_network(const int &n_cluster, const int &R_flag, const int &avoid_resistance_flag, const int &vtk_flag, const vector<double> &voltages, const map<long int, long int> &LMM_cnts, const map<long int, long int> &LMM_gnps, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, Hoshen_Kopelman *HoKo)
{
    //Vectors to extract the backbone
    vector<vector<double> > currents_cnt, currents_gnp;
    
    //Find the "Zero" current of the circuit
    double zero_current;
    //hout << "Find_zero_current" << endl;
    if (!Find_zero_current(n_cluster, R_flag, voltages, LMM_cnts, LMM_gnps, HoKo, cutoffs, points_cnt, radii, points_gnp, gnps, currents_cnt, currents_gnp, zero_current)) {
        hout<<"Error in Determine_backbone_network when calling Find_zero_current"<<endl;
        return 0;
    }
    
    //If there are CNT clusters and there are CNTs in the cluster n_cluster,
    //find the dead branches of the CNTs
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) {
        
        //Find the CNTs in the backbone and dead branches
        //hout << "Find_backbone_and_fractions_cnts" << endl;
        if (!Find_backbone_and_fractions_cnts(n_cluster, avoid_resistance_flag, vtk_flag, zero_current, currents_cnt, points_cnt, radii, HoKo)) {
            hout<<"Error in Determine_backbone_network when calling Find_backbone_and_fractions_cnts"<<endl;
            return 0;
        }
        
    }
    
    //If there are GNP clusters and there are GNPs in the cluster n_cluster, find the dead GNPs
    if (HoKo->clusters_gnp.size() && HoKo->clusters_gnp[n_cluster].size()) {
        
        //Find the GNPs in the backbone and the dead ones
        //hout << "Find_backbone_and_fractions_gnps" << endl;
        if (!Find_backbone_and_fractions_gnps(n_cluster, avoid_resistance_flag, vtk_flag, zero_current, currents_gnp, structure_gnp, gnps, HoKo)) {
            hout<<"Error in Determine_backbone_network when calling Find_backbone_and_fractions_gnps"<<endl;
            return 0;
        }
    }
    
    //Check if resistance will be calcualted to remove the junctions with dead particles
    if (!avoid_resistance_flag) {
        
        //Remove CNT-CNT junctions with dead particles
        //hout << "Remove_junctions_with_dead_particles_cnts" << endl;
        if (!Remove_junctions_with_dead_particles_cnts(n_cluster, HoKo)) {
            hout<<"Error in Determine_backbone_network when calling Remove_junctions_with_dead_particles_cnts"<<endl;
            return 0;
        }
        
        //Remove GNP-GNP junctions with dead particles
        //hout << "Remove_junctions_with_dead_particles_gnps" << endl;
        if (!Remove_junctions_with_dead_particles_gnps(n_cluster, structure_gnp, gnps, HoKo)) {
            hout<<"Error in Determine_backbone_network when calling Remove_junctions_with_dead_particles_gnps"<<endl;
            return 0;
        }
        
        //Remove CNT-GNP junctions with dead particles
        //hout << "Remove_junctions_with_dead_particles_mixed" << endl;
        if (!Remove_junctions_with_dead_particles_mixed(n_cluster, structure_gnp, gnps, HoKo)) {
            hout<<"Error in Determine_backbone_network when calling Remove_junctions_with_dead_particles_mixed"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function finds the "zero current" of the electric circuit
int Backbone_Network::Find_zero_current(const int &n_cluster, const int &R_flag, const vector<double> &voltages, const map<long int, long int> &LMM_cnts, const map<long int, long int> &LMM_gnps, Hoshen_Kopelman *HoKo, const Cutoff_dist &cutoffs, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, vector<vector<double> > &currents_cnt, vector<vector<double> > &currents_gnp, double &zero_current)
{
    //Variable to store the maximum current
    double I_max = Zero;
    
    //Vector to store all the currents
    //vector<double> currents;
    
    //Check if there is any CNT cluster and the current n_cluster has CNTs
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) {
        
        //Initialize the vector of currents for the CNTs
        currents_cnt.assign(HoKo->clusters_cnt[n_cluster].size(), vector<double>());
        
        //Iterate over all CNTs in the cluster
        for (int i = 0; i < (int)HoKo->clusters_cnt[n_cluster].size(); i++) {
            
            //Get the current CNT number
            int CNTi = HoKo->clusters_cnt[n_cluster][i];
            
            //Iterator pointing at the first point in the element set for CNTi
            set<long int>::iterator it = HoKo->elements_cnt[CNTi].begin();
            
            //Get the first point of the element
            long int P1 = *it;
            
            //Get the node of the first point
            long int node1 = LMM_cnts.at(P1);
            
            //Iterate over the elements of CNTi
            for (it++ ; it != HoKo->elements_cnt[CNTi].end(); it++) {
                
                //Get the secont point of the element
                long int P2 = *it;
                
                //Get the second node of the element
                long int node2 = LMM_cnts.at(P2);
                
                //Calcualte the current as the difference between voltages
                //If R_flag is added, then resistances need to be calculated
                double I = abs(voltages[node1] - voltages[node2]);
                
                //Add to the vector of currents
                //currents.push_back(I);
                
                //Check if I is the maximum current so far
                if (I > I_max) {
                    I_max = I;
                }
                
                //Add current to the vector of element currents
                currents_cnt[i].push_back(I);
                
                //Update P1 and node1 for the next iteration
                P1 = P2;
                node1 = node2;
            }
        }
    }
    
    //Check if there are any CNT-CNT junctions in the cluster
    if (HoKo->cluster_cnt_junctions.size() && HoKo->cluster_cnt_junctions[n_cluster].size()) {
        
        //Calculate currents in the junctios and update the maximum current if needed
        if (!Zero_current_in_same_particle_junctions_unit_resistor( HoKo->cluster_cnt_junctions[n_cluster], HoKo->junctions_cnt, voltages, LMM_cnts, I_max)) {
        //if (!Zero_current_in_same_particle_junctions_unit_resistor_test( HoKo->cluster_cnt_junctions[n_cluster], HoKo->junctions_cnt, voltages, LMM_cnts, I_max, currents)) {
            hout<<"Error in Find_zero_current when calling Zero_current_in_same_particle_junctions_unit_resistor (CNTs)"<<endl;
        }
    }
    
    //Check if there is any GNP cluster and the current n_cluster has GNPs
    if (HoKo->clusters_gnp.size() && HoKo->clusters_gnp[n_cluster].size()) {
        
        //Initialize the vector of currents for the GNPs
        currents_gnp.assign(HoKo->clusters_gnp[n_cluster].size(), vector<double>());
        
        //Iterate over all CNTs in the cluster
        for (int i = 0; i < (int)HoKo->clusters_gnp[n_cluster].size(); i++) {
            
            //Get current GNP
            int GNPi = HoKo->clusters_gnp[n_cluster][i];
            
            //Set the size of currents_gnp[i]
            currents_gnp[i].assign(gnps[GNPi].triangulation.size(), 0);
            
            //Go trough all triangulation edges of GNPi
            for (int j = 0; j < (int)gnps[GNPi].triangulation.size(); j++) {
                
                //Get the to triangulation vertices
                long int P1 = gnps[GNPi].triangulation[j].v1;
                long int P2 = gnps[GNPi].triangulation[j].v2;
                
                //Get the nodes
                long int node1 = LMM_gnps.at(P1);
                long int node2 = LMM_gnps.at(P2);
                
                //Calcualte the current as the difference between voltages
                //If R_flag is added, then resistances need to be calculated
                double I = abs(voltages[node1] - voltages[node2]);
                
                //Add to the vector of currents
                //currents.push_back(I);
                
                //Check if I is the maximum current so far
                if (I > I_max) {
                    I_max = I;
                }
                
                //Add current to the vector of element currents
                currents_gnp[i][j] = I;
            }
        }
        
    }
    
    //Check if there are any GNP-GNP junctions in the cluster
    if (HoKo->cluster_gnp_junctions.size() && HoKo->cluster_gnp_junctions[n_cluster].size()) {
        
        //Calculate currents in the junctios and update the maximum current if needed
        if (!Zero_current_in_same_particle_junctions_unit_resistor( HoKo->cluster_gnp_junctions[n_cluster], HoKo->junctions_gnp, voltages, LMM_gnps, I_max)) {
        //if (!Zero_current_in_same_particle_junctions_unit_resistor_test( HoKo->cluster_gnp_junctions[n_cluster], HoKo->junctions_gnp, voltages, LMM_gnps, I_max, currents)) {
            hout<<"Error in Find_zero_current when calling Zero_current_in_same_particle_junctions_unit_resistor (GNPs)"<<endl;
        }
    }
    
    //Check if there are any CNT-GNP junctions in the cluster
    if (HoKo->cluster_mix_junctions.size() && HoKo->cluster_mix_junctions[n_cluster].size()) {
        
        //Iterate over the junctions in the cluster
        for (int i = 0; i < (int)HoKo->cluster_mix_junctions[n_cluster].size(); i++) {
            
            //Get current junction
            int junc = HoKo->cluster_mix_junctions[n_cluster][i];
            
            //Get the points in the junction
            long int P1 = HoKo->junctions_mixed[junc].P1;
            long int P2 = HoKo->junctions_mixed[junc].P2;
            
            //Get the node numbers
            long int node1 = LMM_cnts.at(P1);
            long int node2 = LMM_gnps.at(P2);
            
            //Calcualte the current as the difference between voltages since
            //unit resistors are assumed
            double I = abs(voltages[node1] - voltages[node2]);
            
            //Add to the vector of currents
            //currents.push_back(I);
            
            //Check if I is the maximum current so far
            if (I > I_max) {
                I_max = I;
            }
        }
    }
    
    //Sort currents
    //sort(currents.begin(),currents.end());
    
    //The error cutoff seems to work well with a drop in 9 orders of magnitude of the current. So that is how the cutoff is set.
    //This idea comes from Li and Chou's paper of the DEA in which using a voltage of 1V, a drop in 9 orders of magnitude
    //in the current gave good results.
    //zero_cutoff = currents.back()*1e-9;
    zero_current = I_max*1e-9;
    //hout << "zero_current=" << zero_current << endl;
    
    /* /
    Printer P; string filename;
    if (R_flag) {
        filename = "currents_R_" + to_string(n_cluster) + ".txt";
        P.Print_1d_vec(currents, filename);
    } else {
        filename = "currents_" + to_string(n_cluster) + ".txt";
        P.Print_1d_vec(currents, filename);
    }
    hout << "SAVED CURRENTS FILE " << filename << endl;//*/
    
    return 1;
}
//This function finds the maximum current within the junctions of the same type, i.e.,
//CNT-CNT and GNP-GNP junctions, when unit resitors are used
int Backbone_Network::Zero_current_in_same_particle_junctions_unit_resistor(const vector<int> &junction_cluster, const vector<Junction> &junctions, const vector<double> &voltages, const map<long int, long int> &LMM, double &I_max)
{
    //Iterate over the junctions in the cluster
    for (int i = 0; i < (int)junction_cluster.size(); i++) {
        
        //Get current junction
        int junc = junction_cluster[i];
        
        //Get the points in the junction
        long int P1 = junctions[junc].P1;
        long int P2 = junctions[junc].P2;
        
        //Get the node numbers
        long int node1 = LMM.at(P1);
        long int node2 = LMM.at(P2);
        
        //Calcualte the current as the difference between voltages since
        //unit resistors are assumed
        double I = abs(voltages[node1] - voltages[node2]);
        
        //Check if I is the maximum current so far
        if (I > I_max) {
            I_max = I;
        }
    }
    
    return 1;
}
//Same function as the one above, but this adds the calculated current to the vector of currents
int Backbone_Network::Zero_current_in_same_particle_junctions_unit_resistor_test(const vector<int> &junction_cluster, const vector<Junction> &junctions, const vector<double> &voltages, const map<long int, long int> &LMM, double &I_max, vector<double> &currents)
{
    //Iterate over the junctions in the cluster
    for (int i = 0; i < (int)junction_cluster.size(); i++) {
        
        //Get current junction
        int junc = junction_cluster[i];
        
        //Get the points in the junction
        long int P1 = junctions[junc].P1;
        long int P2 = junctions[junc].P2;
        
        //Get the node numbers
        long int node1 = LMM.at(P1);
        long int node2 = LMM.at(P2);
        
        //Calcualte the current as the difference between voltages since
        //unit resistors are assumed
        double I = abs(voltages[node1] - voltages[node2]);
        
        //Add to the vector of currents
        currents.push_back(I);
        
        //Check if I is the maximum current so far
        if (I > I_max) {
            I_max = I;
        }
    }
    
    return 1;
}
//This function finds the backbone and calculates the fractions of percolated families
//If needed, VTK files are generated and HoKo is updated for calculating electrical resistance
int Backbone_Network::Find_backbone_and_fractions_cnts(const int &n_cluster, const int &avoid_resistance_flag, const int &vtk_flag, const double &zero_current, const vector<vector<double> > &currents_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, Hoshen_Kopelman *HoKo)
{
    //Vectors to export vtk visualization files
    vector<vector<long int> > dead_branches_idx, backbone_idx;
    
    //Initialize the vectors for the vtk visualization files if they are needed
    if (vtk_flag) {
        dead_branches_idx.assign(HoKo->clusters_cnt[n_cluster].size(), vector<long int>());
        backbone_idx.assign(HoKo->clusters_cnt[n_cluster].size(), vector<long int>());
    }
    
    //Iterate over the CNTs in the cluster
    //Iterate from the last element to the first in order to be able to delete elements
    for (int i = (int)HoKo->clusters_cnt[n_cluster].size()-1; i >= 0 ; i--) {
        
        //Get the current CNT
        //hout << "n_cluster=" << n_cluster << " HoKo->clusters_cnt.size=" << HoKo->clusters_cnt.size() << endl;
        int CNTi = HoKo->clusters_cnt[n_cluster][i];
        //hout << "CNTi=" << CNTi << endl;
        
        //Indices that determine the segment of the CNT that belongs to the backbone
        int idx1 = -1;
        int idx2 = -1;
        
        //Iterate over the currents vector
        for (int j = 0; j < (int)currents_cnt[i].size(); j++) {
            
            //Check if the current is above  zero current
            //hout << "j=" << j << endl;
            if (currents_cnt[i][j] - zero_current > Zero) {
                
                //This CNT segment is part of the backbone
                //Save the first index if this is the first backbone segment found
                //This happens when idx is still -1
                if (idx1 == -1) {
                    idx1 = j;
                }
            }
            else {
                
                //This CNT segment is part of the dead branches
                //Save the second index, if this is the first dead segment after
                //a backbone segment
                //This happens if idx1 is not -1 and idx2 is still -1
                if (idx1 != -1 && idx2 == -1) {
                    idx2 = j;
                    //Terminate the loop, there is no need to continue looking for
                    //dead branches after a percolated segment was found
                    break;
                }
            }
        }
        
        //Calcualte the volumes of the CNT segments for the dead branches and backbone, then
        //add them to the class variables
        //If needed, add indices to the vectors needed to export VTK files and remove points
        //from the elements vector to calculate electrical resistance in a further step
        //hout << "Calculate_cnt_volumes" << endl;
        if (!Calculate_cnt_volumes(n_cluster, avoid_resistance_flag, vtk_flag, CNTi, idx1, idx2, i, points_cnt, radii[CNTi], HoKo, dead_branches_idx, backbone_idx)) {
            hout<<"Error in Find_backbone_and_fractions_cnts when calling Calculate_cnt_volumes"<<endl;
            return 0;
        }
        
        //Check if the whole CNT is dead and remove it from the cluster in case
        //electrical resistance is to be calculated
        //hout << "Check if the whole CNT is dead" << endl;
        if (!avoid_resistance_flag && HoKo->elements_cnt[CNTi].empty()) {
            HoKo->clusters_cnt[n_cluster].erase(HoKo->clusters_cnt[n_cluster].begin()+i);
        }
    }
    
    //Export the visualization files if needed
    if (vtk_flag) {
        
        //VTK export object
        VTK_Export VTK_E;
        
        //Generate filenames
        string str_bb = "backbone_C" + to_string(n_cluster) + "_F" + to_string(HoKo->family[n_cluster]) + "_cnts.vtk";
        string str_db = "dead_branches_C" + to_string(n_cluster) + "_F" + to_string(HoKo->family[n_cluster]) + "_cnts.vtk";
        
        //Export CNTs in the backbone
        if (!VTK_E.Export_from_cnt_indices(points_cnt, backbone_idx, str_bb)) {
            hout<<"Error in Export_vtk_files_for_backbone_cnts when calling VTK_E.Export_from_cnt_indices (backbone)"<<endl;
            return 0;
        }
        
        //Export dead branches
        if (!VTK_E.Export_from_cnt_indices(points_cnt, dead_branches_idx, str_db)) {
            hout<<"Error in Export_vtk_files_for_backbone_cnts when calling VTK_E.Export_from_cnt_indices (dead branches)"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function calcualtes the volumes of the segments of a CNT that belong to dead branches
//and backbone
int Backbone_Network::Calculate_cnt_volumes(const int &n_cluster, const int &avoid_resistance_flag, const int &vtk_flag, const int &CNTi, const int &idx1, const int &idx2, const int& branch_idx, const vector<Point_3D> &points_cnt, const double &radius, Hoshen_Kopelman *HoKo, vector<vector<long int> >& dead_branches_idx, vector<vector<long int> >& backbone_idx)
{
    //Get the two points at the beginning and end of the elements vector
    //hout << "HoKo->elements_cnt.size=" << HoKo->elements_cnt.size() << endl;
    long int P_start = *(HoKo->elements_cnt[CNTi].begin());
    long int P_end = *(HoKo->elements_cnt[CNTi].rbegin());
    //hout << "P_start=" << P_start << " P_end=" << P_end << endl;
    
    //Variables to store the volumes for dead branches and backbone
    double volume_backbone = 0, volume_dead = 0;
    
    //Get the family number
    int family = HoKo->family[n_cluster];
    //hout << "family=" << family << endl;
    
    //If idx1 is still -1, then all the CNT is dead
    //If idx1 is not -1, then
    //  If idx2 is -1, then starting at idx1 the rest of the CNT is part of the backbone
    //  If idx2 is not -1, then only the segment between idx1 and idx2 is part of the backbone
    if (idx1 == -1) {
        
        //Calculate the volume of the whole CNT
        //hout << "CNT_volume_between_two_points" << endl;
        if (!CNT_volume_between_two_points(P_start, P_end, radius, points_cnt, volume_dead)) {
            hout<<"Error in Calculate_cnt_volumes when calling CNT_volume_between_two_points"<<endl;
            return 0;
        }
        
        //Add the volume of the whole CNT to the dead branches of the family the
        //cluster belongs to
        dead_branches[family] = dead_branches[family] + volume_dead;
        //hout << "dead_branches[fam="<<family<<"]=" << dead_branches[family] << endl;
        
        //If visualization files are needed, add the indices to the vector
        if (vtk_flag) {
            dead_branches_idx[branch_idx].push_back(P_start);
            dead_branches_idx[branch_idx].push_back(P_end);
        }
        
        //Remove all points in the elements if resistance will be calcualted
        if (!avoid_resistance_flag) {
            HoKo->elements_cnt[CNTi].clear();
        }
    }
    else {
        
        //Get the point at idx1
        set<long int>::iterator it1 = HoKo->elements_cnt[CNTi].begin();
        advance(it1, idx1);
        long int P1 = *it1;
        
        //Get the volumes from P_start to P1, which corresponds to the dead brach
        //hout << "CNT_volume_between_two_points P1=" << P1 << endl;
        if (!CNT_volume_between_two_points(P_start, P1, radius, points_cnt, volume_dead)) {
            hout<<"Error in Calculate_cnt_volumes when calling CNT_volume_between_two_points"<<endl;
            return 0;
        }
        
        //Add the volume of the the dead branch to the corresponding family
        dead_branches[family] = dead_branches[family] + volume_dead;
        //hout << "dead_branches[fam=" << family << "]=" << dead_branches[family] << endl;
        
        //If visualization files are needed, add the indices to the vector
        if (vtk_flag) {
            dead_branches_idx[branch_idx].push_back(P_start);
            dead_branches_idx[branch_idx].push_back(P1);
        }
        
        //Check the value of idx2
        if (idx2 == -1) {
            
            //Get the volumes from P1 to P_end, which corresponds to the backbone
            //hout << "CNT_volume_between_two_points idx2=-1" << endl;
            if (!CNT_volume_between_two_points(P1, P_end, radius, points_cnt, volume_backbone)) {
                hout<<"Error in Calculate_cnt_volumes when calling CNT_volume_between_two_points"<<endl;
                return 0;
            }
            
            //If visualization files are needed, add the indices to the vector
            if (vtk_flag) {
                backbone_idx[branch_idx].push_back(P1);
                backbone_idx[branch_idx].push_back(P_end);
            }
        }
        else {
            
            //There are two dead branches, so get the point at idx2 for the second dead branch
            it1 = HoKo->elements_cnt[CNTi].begin();
            advance(it1, idx2);
            long int P2 = *it1;
            
            //Get the volumes from P1 to P2, which corresponds to the backbone
            //hout << "CNT_volume_between_two_points P2=" << P2 << endl;
            if (!CNT_volume_between_two_points(P1, P2, radius, points_cnt, volume_backbone)) {
                hout<<"Error in Calculate_cnt_volumes when calling CNT_volume_between_two_points"<<endl;
                return 0;
            }
            
            //If visualization files are needed, add the indices to the vector
            if (vtk_flag) {
                backbone_idx[branch_idx].push_back(P1);
                backbone_idx[branch_idx].push_back(P2);
            }
            
            //Get the volumes from P2 to P_end, which corresponds to the second dead brach
            //hout << "CNT_volume_between_two_points" << endl;
            if (!CNT_volume_between_two_points(P2, P_end, radius, points_cnt, volume_dead)) {
                hout<<"Error in Calculate_cnt_volumes when calling CNT_volume_between_two_points"<<endl;
                return 0;
            }
            
            //Add the volume of the the dead branch to the corresponding family
            dead_branches[family] = dead_branches[family] + volume_dead;
            //hout << "dead_branches[fam=" << family << "]=" << dead_branches[family] << endl;
            
            //Erase the points for dead branches if the resistance will be calculated
            if (!avoid_resistance_flag) {
                //Remove the points after P2 (idx2) from the elements vector
                //it1 already points to P2, so increase it by one to point to the next point
                it1++;
                HoKo->elements_cnt[CNTi].erase(it1, HoKo->elements_cnt[CNTi].end());
            }
            
            //If visualization files are needed, add the indices to the vector
            if (vtk_flag) {
                dead_branches_idx[branch_idx].push_back(P2);
                dead_branches_idx[branch_idx].push_back(P_end);
            }
        }
        
        //Add the volume to the percolated family
        volumes_cnt[family] = volumes_cnt[family] + volume_backbone;
        //hout << "volumes_cnt[fam=" << family << "]=" << volumes_cnt[family] << endl;
        
        //Erase the points for dead branches if the resistance will be calculated
        //Remove all points before P1 (idx1) from the elements vector, except when idx1 is zero
        if (!avoid_resistance_flag && idx1 > 0) {
            it1 = HoKo->elements_cnt[CNTi].begin();
            advance(it1, idx1);
            HoKo->elements_cnt[CNTi].erase(HoKo->elements_cnt[CNTi].begin(), it1);
        }
    }

    return 1;
}
//This function calculate the volume of a CNT between two points
int Backbone_Network::CNT_volume_between_two_points(const long int &P1, const long int P2, const double &radius, const vector<Point_3D> &points_cnt, double &volume)
{
    //Variable to store the length of the segment between P1 and P2
    double length = 0;
    
    //Iterate from P1 up to P2
    for (long int i = P1; i < P2; i++) {
        
        //Add the distance from i to i+1 to the total length of the segment
        length = length + points_cnt[i].distance_to(points_cnt[i+1]);
    }
    
    //Multiply the length by the cross-sectional area to obtain the volume
    volume = length*PI*radius*radius;
    
    return 1;
}
//This function finds the backbone and calculates the fractions of percolated families of GNPs
//If needed, VTK files are generated and HoKo is updated for calculating electrical resistance
int Backbone_Network::Find_backbone_and_fractions_gnps(const int &n_cluster, const int &avoid_resistance_flag, const int &vtk_flag, const double &zero_current, const vector<vector<double> > &currents_gnp, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, Hoshen_Kopelman *HoKo)
{
    //Vectors to export vtk visualization files
    vector<int> dead_gnps_idx, backbone_gnps_idx;
    
    //Get the family
    //hout << "n_cluster=" << n_cluster << " HoKo->family.size=" << HoKo->family.size() << endl;
    int family = HoKo->family[n_cluster];
    //hout << "family=" << family << endl;
    
    //Iterate over the GNPs in the cluster
    //Iterate from the last GNP to the first in order to be able to delete them
    for (int i = (int)HoKo->clusters_gnp[n_cluster].size()-1; i >= 0 ; i--) {
        
        //hout << "i=" << i << endl;
        //Get current GNP
        int GNPi = HoKo->clusters_gnp[n_cluster][i];
        
        //non-zero current flag
        bool nz_flag = false;
        
        //Iterate over the currents of the triangulation edges, which are saved in
        //currents_gnp[i]
        for (int j = 0; j < (int)currents_gnp[i].size(); j++) {
            
            //hout << "j=" << j << endl;
            //Check if current in triangulation edge is non-zero
            if (currents_gnp[i][j] > zero_current) {
                
                //hout << "currents_gnp[i][j]=" << currents_gnp[i][j] << endl;
                //There is current flowign through GNPi, so it is part of the backbone
                //Set the non-zero current flag to true
                nz_flag = true;
                
                //Break the loop as it is not necessary to check for more non-zero
                //currents to conclude GNPi is part of the backbone
                break;
            }
        }
        
        //Check if the volume of the GNP should be added to the percolated or dead part
        //hout << "nz_flag=" << endl;
        if (nz_flag) {
            
            //Add GNP volume to percolated part
            volumes_gnp[family] = volumes_gnp[family] + gnps[GNPi].volume;
            
            //If visualization files are needed, add the GNP to backbone
            if (vtk_flag) {
                backbone_gnps_idx.push_back(GNPi);
            }
        }
        else {
            
            //Add GNP volume to dead part
            dead_gnps[family] = dead_gnps[family] + gnps[GNPi].volume;
            
            //If visualization files are needed, add the GNP to the dead ones
            if (vtk_flag) {
                dead_gnps_idx.push_back(GNPi);
            }
        }
        
        //If the electrical resistance will be calcualted, remove all triangulation edges of GNPi
        //hout << "avoid_resistance_flag" << endl;
        if (!avoid_resistance_flag) {
            
            //Clear the triangulation
            gnps[GNPi].triangulation.clear();
            
            //Check if current GNP is dead and, if so, remove it form the cluster and
            //remove all points from the structure
            if (!nz_flag) {
                HoKo->clusters_gnp[n_cluster].erase(HoKo->clusters_gnp[n_cluster].begin()+i);
                structure_gnp[GNPi].clear();
            }
        }
    }
    
    //Check if visualization files are needed
    if (vtk_flag) {
        //VTK export object
        VTK_Export VTK_E;
        
        //Generate filenames for CNT files
        string str_bb = "backbone_C" + to_string(n_cluster) + "_F" + to_string(HoKo->family[n_cluster]) + "_gnps.vtk";
        string str_dead = "dead_C" + to_string(n_cluster) + "_F" + to_string(HoKo->family[n_cluster]) + "_gnps.vtk";
        
        //Export GNPs in the backbone
        if (!VTK_E.Export_gnps_in_cluster(gnps, backbone_gnps_idx, str_bb)) {
            hout<<"Error in Find_backbone_and_fractions_gnps when calling VTK_E.Export_gnps_in_cluster"<<endl;
            return 0;
        }
        
        //Export dead GNPs
        if (!VTK_E.Export_gnps_in_cluster(gnps, dead_gnps_idx, str_dead)) {
            hout<<"Error in Find_backbone_and_fractions_gnps when calling VTK_E.Export_gnps_in_cluster"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function removes the junctions with dead particles for CNT-CNT junctions
int Backbone_Network::Remove_junctions_with_dead_particles_cnts(const int &n_cluster, Hoshen_Kopelman *HoKo)
{
    //Check if there are junctions in the cluster
    if (HoKo->cluster_cnt_junctions.size() && HoKo->cluster_cnt_junctions[n_cluster].size()) {
        
        //Iterate over the junctions in the cluster
        //Iterate from the last element to be able delete elements
        for (int i = (int)HoKo->cluster_cnt_junctions[n_cluster].size() - 1; i >= 0 ; i--) {
            
            //Get the current junction
            int junc = HoKo->cluster_cnt_junctions[n_cluster][i];
            
            //Get the CNT numbers
            int CNT1 = HoKo->junctions_cnt[junc].N1;
            int CNT2 = HoKo->junctions_cnt[junc].N2;
            
            //Check if CNTs are dead
            if (HoKo->elements_cnt[CNT1].empty() || HoKo->elements_cnt[CNT2].empty()) {
                
                //At least one of the CNTs is dead, so remove the junction
                HoKo->cluster_cnt_junctions[n_cluster].erase(HoKo->cluster_cnt_junctions[n_cluster].begin()+i);
            }
            else {
                
                //Both CNTs are in the backbone
                //Get the points of the junction
                long int P1 = HoKo->junctions_cnt[junc].P1;
                long int P2 = HoKo->junctions_cnt[junc].P2;
                
                //Check if both points are in the elements vector
                if (HoKo->elements_cnt[CNT1].find(P1) == HoKo->elements_cnt[CNT1].end() ||
                    HoKo->elements_cnt[CNT2].find(P2) == HoKo->elements_cnt[CNT2].end()) {
                    
                    //At least one of the two points is not in the backbone
                    //Then, remove this junction from the cluster of junctions
                    HoKo->cluster_cnt_junctions[n_cluster].erase(HoKo->cluster_cnt_junctions[n_cluster].begin()+i);
                }
            }
        }
    }
    
    return 1;
}
//This function removes the junctions with dead particles for GNP-GNP junctions
int Backbone_Network::Remove_junctions_with_dead_particles_gnps(const int &n_cluster, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, Hoshen_Kopelman *HoKo)
{
    //Check if there are junctions in the cluster
    if (HoKo->cluster_gnp_junctions.size() && HoKo->cluster_gnp_junctions[n_cluster].size()) {
        
        //Iterate over the junctions in the cluster
        //Iterate from the last element to be able delete elements
        for (int i = (int)HoKo->cluster_gnp_junctions[n_cluster].size() - 1; i >= 0 ; i--) {
            
            //Get the current junction
            int junc = HoKo->cluster_gnp_junctions[n_cluster][i];
            
            //Get the GNP numbers
            int GNP1 = HoKo->junctions_gnp[junc].N1;
            int GNP2 = HoKo->junctions_gnp[junc].N2;
            
            //Check if there are points in the structure vector
            if (structure_gnp[GNP1].empty() || structure_gnp[GNP2].empty()) {
                
                //At least one GNP is dead, so remove the junction from the cluster
                HoKo->cluster_gnp_junctions[n_cluster].erase(HoKo->cluster_gnp_junctions[n_cluster].begin()+i);
                
                //If GNP1 is not dead, remove point 1 of the junction from the structure
                if (!structure_gnp[GNP1].empty()) {
                    
                    //Get the point in GNP1
                    long int P1 = HoKo->junctions_gnp[junc].P1;
                    if (!Remove_point_from_vector(P1, structure_gnp[GNP1])) {
                        hout<<"Error in Remove_junctions_with_dead_particles_gnps when calling Remove_point_from_vector (1)"<<endl;
                        return 0;
                    }
                }
                
                //If GNP2 is not dead, remove point 2 of the junction from the structure
                if (!structure_gnp[GNP2].empty()) {
                    
                    //Get the point in GNP2
                    long int P2 = HoKo->junctions_gnp[junc].P2;
                    if (!Remove_point_from_vector(P2, structure_gnp[GNP2])) {
                        hout<<"Error in Remove_junctions_with_dead_particles_gnps when calling Remove_point_from_vector (2)"<<endl;
                        return 0;
                    }
                }
            }
        }
    }
    
    return 1;
}
//This function removes a point from a vector
int Backbone_Network::Remove_point_from_vector(const long int &P, vector<long int> &structure_i)
{
    for (int j = 0; j < (int)structure_i.size(); j++) {
        if (structure_i[j] == P) {
            structure_i.erase(structure_i.begin()+j);
            break;
        }
    }
    
    return 1;
}
//This function removes the junctions with dead particles for CNT-GNP junctions
int Backbone_Network::Remove_junctions_with_dead_particles_mixed(const int &n_cluster, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, Hoshen_Kopelman *HoKo)
{
    //Check if there are junctions in the cluster
    if (HoKo->cluster_mix_junctions.size() && HoKo->cluster_mix_junctions[n_cluster].size()) {
        
        //Iterate over the junctions in the cluster
        //Iterate from the last element to be able delete elements
        for (int i = (int)HoKo->cluster_mix_junctions[n_cluster].size() - 1; i >= 0 ; i--) {
            
            //Get the current junction
            int junc = HoKo->cluster_mix_junctions[n_cluster][i];
            
            //Get the particle numbers
            int CNT1 = HoKo->junctions_cnt[junc].N1;
            int GNP2 = HoKo->junctions_cnt[junc].N2;
            
            //Check if CNT or GNP is dead
            if (HoKo->elements_cnt[CNT1].empty() || structure_gnp[GNP2].empty()) {
                
                //At least one of the particles is dead, so remove the junction
                HoKo->cluster_mix_junctions[n_cluster].erase(HoKo->cluster_mix_junctions[n_cluster].begin()+i);
                
                //If GNP2 is not dead, remove point 2 of the junction from the structure
                if (!structure_gnp[GNP2].empty()) {
                    
                    //Get the point in GNP2 and remove it from the structure
                    long int P2 = HoKo->junctions_mixed[junc].P2;
                    if (!Remove_point_from_vector(P2, structure_gnp[GNP2])) {
                        hout<<"Error in Remove_junctions_with_dead_particles_gnps when calling Remove_point_from_vector (2)"<<endl;
                        return 0;
                    }
                }
            }
            else {
                
                //Both particles are in the backbone
                //Get the CNT point of the junction
                long int P1 = HoKo->junctions_mixed[junc].P1;
                
                //Check if P1 is in the elements vector
                if (HoKo->elements_cnt[CNT1].find(P1) == HoKo->elements_cnt[CNT1].end()) {
                    
                    //Although both particles are in the backbone, P1 is not
                    //Then, remove this junction from the cluster of junctions
                    HoKo->cluster_mix_junctions[n_cluster].erase(HoKo->cluster_mix_junctions[n_cluster].begin()+i);
                    
                    //Get the point in GNP2 and remove it from the structure
                    long int P2 = HoKo->junctions_mixed[junc].P2;
                    if (!Remove_point_from_vector(P2, structure_gnp[GNP2])) {
                        hout<<"Error in Remove_junctions_with_dead_particles_gnps when calling Remove_point_from_vector (2)"<<endl;
                        return 0;
                    }
                }
            }
        }
    }
    
    return 1;
}
