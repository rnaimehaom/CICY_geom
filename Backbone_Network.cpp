//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Determine the backbone network and dead branches in the percolation network
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Backbone_Network.h"


int Backbone_Network::Determine_backbone_network(const int &family, const int &n_cluster, const int &R_flag, const int &avoid_resistance_flag, const int &vtk_flag, const vector<double> &voltages, const map<long int, long int> &LMM_cnts, const map<long int, long int> &LMM_gnps, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, Hoshen_Kopelman *HoKo)
{
    //Vectors to extract the backbone
    vector<vector<double> > currents_cnt, currents_gnp;
    
    //Find the "Zero" current of the circuit
    double zero_current;
    if (!Find_zero_current(n_cluster, R_flag, voltages, LMM_cnts, LMM_gnps, HoKo, cutoffs, points_cnt, radii, points_gnp, gnps, currents_cnt, currents_gnp, zero_current)) {
        hout<<"Error in Determine_backbone_network when calling Find_zero_current"<<endl;
        return 0;
    }
    
    //If there are CNT clusters and there are CNTs in the cluster n_cluster,
    //find the dead branches of the CNTs
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) {
        
        //Find the CNTs in the backbone and dead branches
        if (Find_backbone_and_fractions_cnts(family, n_cluster, avoid_resistance_flag, vtk_flag, zero_current, currents_cnt, points_cnt, radii, HoKo)) {
            hout<<"Error in Determine_backbone_network when calling Find_backbone_and_fractions_cnts"<<endl;
            return 0;
        }
        
    }
    
    //If there are GNP clusters and there are GNPs in the cluster n_cluster, find the dead GNPs
    if (HoKo->clusters_gnp.size() && HoKo->clusters_gnp[n_cluster].size()) {
        
        //Find the GNPs in the backbone and the dead ones
        if (!Find_backbone_and_fractions_gnps(family, n_cluster, avoid_resistance_flag, vtk_flag, zero_current, currents_gnp, structure_gnp, gnps, HoKo)) {
            hout<<"Error in Determine_backbone_network when calling Find_backbone_and_fractions_gnps"<<endl;
            return 0;
        }
    }
    
    //Check if resistance will be calcualted to remove the junctions with dead particles
    if (!avoid_resistance_flag) {
        
        //Remove CNT-CNT junctions with dead particles
        if (!Remove_junctions_with_dead_particles_cnts(n_cluster, HoKo)) {
            hout<<"Error in Determine_backbone_network when calling Remove_junctions_with_dead_particles_cnts"<<endl;
            return 0;
        }
        
        //Remove GNP-GNP junctions with dead particles
        if (!Remove_junctions_with_dead_particles_gnps(n_cluster, structure_gnp, gnps, HoKo)) {
            hout<<"Error in Determine_backbone_network when calling Remove_junctions_with_dead_particles_gnps"<<endl;
            return 0;
        }
        
        //Remove CNT-GNP junctions with dead particles
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
    
    /*/
    Printer *P = new Printer;
    if (R_flag) {
        P->Print_1d_vec(currents, "currents_R.txt");
        P->Print_2d_vec(currents_gnp, "currents_gnp_R.txt");
    } else {
        P->Print_1d_vec(currents, "currents.txt");
        P->Print_2d_vec(currents_gnp, "currents_gnp.txt");
    }
    delete P;//*/
    
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
//This function finds the backbone and calculates the fractions of percolated families
//If needed, VTK files are generated and HoKo is updated for calculating electrical resistance
int Backbone_Network::Find_backbone_and_fractions_cnts(const int &family, const int &n_cluster, const int &avoid_resistance_flag, const int &vtk_flag, const double &zero_current, const vector<vector<double> > &currents_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, Hoshen_Kopelman *HoKo)
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
        int CNTi = HoKo->clusters_cnt[n_cluster][i];
        
        //Indices that determine the segment of the CNT that belongs to the backbone
        int idx1 = -1;
        int idx2 = -1;
        
        //Iterate over the currents vector
        for (int j = 0; j < (int)currents_cnt[i].size(); j++) {
            
            //Check if the current is above  zero current
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
        
        //Calcualte the volumes of the CNt segments for the ded branches and backbone, then
        //add them to the class variables
        //If needed, add indices to the vectors needed to export VTK files and remove points
        //from the elements vector to calculate electrical resistance in a further step
        if (!Calculate_cnt_volumes(family, n_cluster, avoid_resistance_flag, vtk_flag, CNTi, idx1, idx2, points_cnt, radii[CNTi], HoKo, dead_branches_idx[i], backbone_idx[i])) {
            hout<<"Error in Find_backbone_and_fractions_cnts when calling Calculate_cnt_volumes"<<endl;
            return 0;
        }
        
        //Check if the whole CNT is dead and remove it from the cluster in case
        //electrical resistance is to be calculated
        if (!avoid_resistance_flag && HoKo->elements_cnt[CNTi].empty()) {
            HoKo->clusters_cnt[n_cluster].erase(HoKo->clusters_cnt[n_cluster].begin()+i);
        }
    }
    
    //Export the visualization files if needed
    if (vtk_flag) {
        //
    }
    
    return 1;
}
//This function calcualtes the volumes of the segments of a CNT that belong to dead branches
//and backbone
int Backbone_Network::Calculate_cnt_volumes(const int &family, const int &n_cluster, const int &avoid_resistance_flag, const int &vtk_flag, const int &CNTi, const int &idx1, const int &idx2, const vector<Point_3D> &points_cnt, const double &radius, Hoshen_Kopelman *HoKo, vector<long int> &dead_branches_i, vector<long int> &backbone_i)
{
    //Get the two points at the beginning and end of the elements vector
    long int P_start = *(HoKo->elements_cnt[CNTi].begin());
    long int P_end = *(HoKo->elements_cnt[CNTi].rbegin());
    
    //Variables to store the volumes for dead branches and backbone
    double volume_backbone = 0, volume_dead = 0;
    
    //If idx1 is still -1, then all the CNT is dead
    //If idx1 is not -1, then
    //  If idx2 is -1, then starting at idx1 the rest of the CNT is part of the backbone
    //  If idx2 is not -1, then only the segment between idx1 and idx2 is part of the backbone
    if (idx1 == -1) {
        
        //Calculate the volume of the whole CNT
        if (!CNT_volume_between_two_points(P_start, P_end, radius, points_cnt, volume_dead)) {
            hout<<"Error in Calculate_cnt_volumes when calling CNT_volume_between_two_points"<<endl;
            return 0;
        }
        
        //Add the volume of the whole CNT to the dead branches of the family the
        //cluster belongs to
        dead_branches[family] = dead_branches[family] + volume_dead;
        
        //If visualization files are needed, add the indices to the vector
        if (vtk_flag) {
            dead_branches_i.push_back(P_start);
            dead_branches_i.push_back(P_end);
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
        if (!CNT_volume_between_two_points(P_start, P1, radius, points_cnt, volume_dead)) {
            hout<<"Error in Calculate_cnt_volumes when calling CNT_volume_between_two_points"<<endl;
            return 0;
        }
        
        //Add the volume of the the dead branch to the corresponding family
        dead_branches[family] = dead_branches[family] + volume_dead;
        
        //If visualization files are needed, add the indices to the vector
        if (vtk_flag) {
            dead_branches_i.push_back(P_start);
            dead_branches_i.push_back(P1);
        }
        
        //Check the value of idx2
        if (idx2 == -1) {
            
            //Get the volumes from P1 to P_end, which corresponds to the backbone
            if (!CNT_volume_between_two_points(P1, P_end, radius, points_cnt, volume_backbone)) {
                hout<<"Error in Calculate_cnt_volumes when calling CNT_volume_between_two_points"<<endl;
                return 0;
            }
            
            //If visualization files are needed, add the indices to the vector
            if (vtk_flag) {
                backbone_i.push_back(P1);
                backbone_i.push_back(P_end);
            }
        }
        else {
            
            //There are two dead branches, so get the point at idx2 for the second dead branch
            it1 = HoKo->elements_cnt[CNTi].begin();
            advance(it1, idx2);
            long int P2 = *it1;
            
            //Get the volumes from P1 to P2, which corresponds to the backbone
            if (!CNT_volume_between_two_points(P1, P2, radius, points_cnt, volume_backbone)) {
                hout<<"Error in Calculate_cnt_volumes when calling CNT_volume_between_two_points"<<endl;
                return 0;
            }
            
            //If visualization files are needed, add the indices to the vector
            if (vtk_flag) {
                backbone_i.push_back(P1);
                backbone_i.push_back(P2);
            }
            
            //Get the volumes from P2 to P_end, which corresponds to the second dead brach
            if (!CNT_volume_between_two_points(P2, P_end, radius, points_cnt, volume_dead)) {
                hout<<"Error in Calculate_cnt_volumes when calling CNT_volume_between_two_points"<<endl;
                return 0;
            }
            
            //Add the volume of the the dead branch to the corresponding family
            dead_branches[family] = dead_branches[family] + volume_dead;
            
            //Erase the points for dead branches if the resistance will be calculated
            if (!avoid_resistance_flag) {
                //Remove the points after P2 (idx2) from the elements vector
                //it1 already points to P2, so increase it by one to point to the next point
                it1++;
                HoKo->elements_cnt[CNTi].erase(it1, HoKo->elements_cnt[CNTi].end());
            }
            
            //If visualization files are needed, add the indices to the vector
            if (vtk_flag) {
                dead_branches_i.push_back(P2);
                dead_branches_i.push_back(P_end);
            }
        }
        
        //Add the volume to the percolated family
        volumes_cnt[family] = volumes_cnt[family] + volume_backbone;
        
        //Erase the points for dead branches if the resistance will be calculated
        //Remove all points before P1 (idx1) from the elements vector, except when idx1 is zero
        if (!avoid_resistance_flag && idx1 > 0) {
            it1 = HoKo->elements_cnt[CNTi].begin();
            advance(it1, idx1-1);
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
//This function generates two VTK files for a cluster of CNTs: one for dead branches
//and one for the backbone
int Export_percolated_and_non_percoalted_clusters()
{
    return 1;
}
//This function finds the backbone and calculates the fractions of percolated families of GNPs
//If needed, VTK files are generated and HoKo is updated for calculating electrical resistance
int Backbone_Network::Find_backbone_and_fractions_gnps(const int &family, const int &n_cluster, const int &avoid_resistance_flag, const int &vtk_flag, const double &zero_current, const vector<vector<double> > &currents_gnp, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps, Hoshen_Kopelman *HoKo)
{
    //Vectors to export vtk visualization files
    vector<int> dead_gnps_idx, backbone_gnps_idx;
    
    //Iterate over the GNPs in the cluster
    //Iterate from the last GNP to the first in order to be able to delete them
    for (int i = (int)HoKo->clusters_gnp[n_cluster].size()-1; i >= 0 ; i--) {
        
        //Get current GNP
        int GNPi = HoKo->clusters_gnp[n_cluster][i];
        
        //non-zero current flag
        bool nz_flag = false;
        
        //Iterate over the currents of the triangulation edges, which are saved in
        //currents_gnp[i]
        for (int j = 0; j < (int)currents_gnp[i].size(); j++) {
            
            //Check if current in triangulation edge is non-zero
            if (currents_gnp[i][j] > zero_current) {
                
                //There is current flowign through GNPi, so it is part of the backbone
                //Set the non-zero current flag to true
                nz_flag = true;
                
                //Break the loop as it is not necessary to check for more non-zero
                //currents to conclude GNPi is part of the backbone
                break;
            }
        }
        
        //Check if the volume of the GNP should be added to the percolated or dead part
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
        
        //Export all GNPs in the backbone and the dead ones
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
//
//Deprecated:
//
int Backbone_Network::Determine_backbone_network(const int &family, const int &n_cluster, const int &R_flag, Direct_Electrifying *DEA, Hoshen_Kopelman *HoKo, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_in_gnp, const vector<GCH> &hybrid_particles, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_percolated_indices, vector<vector<int> > &all_dead_gnps, vector<vector<int> > &all_percolated_gnps)
{
    //First check if the CNT cluster has only one CNT, in that case there are no dead branches and the backbone extraction can be simplified a lot
    //Otherwise look for dead branches
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size() == 1 && !HoKo->clusters_gch.size()) {
        
        //When there is only one percolated CNT, then the percolated_indices only contains
        //the first and last point of the CNT and dead_indices remains empty
        int CNT = HoKo->clusters_cnt[n_cluster].front();
        //This varibale is used to initialize the vectors below
        vector<long int> empty;
        percolated_indices.push_back(empty);
        percolated_indices.back().push_back(structure[CNT].front());
        percolated_indices.back().push_back(structure[CNT].back());
        //Initialize dead_indices as it needs to have the same size as percolated_indices
        dead_indices.push_back(empty);
        
    } else {
        
        //Current vectors to determine conducting segments and conducting GNPs
        vector<vector<double> > currents_cnt, currents_gnp;
        
        //Define the 'Zero-cutoff' of the system
        double zero_cutoff = Zero_current(n_cluster, R_flag, DEA, HoKo, electric_param, cutoffs, points_in, points_in_gnp, radii, hybrid_particles, currents_cnt, currents_gnp);
        
        //If there are CNT clusters and there are CNTs in the cluster n_cluster, find the dead branches of the CNTs
        if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) {
            
            //Find the dead branches of each CNT in the cluster. Save the information on the vectors dead_indices and percolated_indices
            if (!Find_dead_branches(zero_cutoff, currents_cnt, HoKo->clusters_cnt[n_cluster], DEA->elements)) {
                hout << "Error in Determine_backbone_network when calling Find_dead_branches" << endl;
                return 0;
            }
            
            //Add the indices of percolated and non percolated CNTs
            if (!Add_indices_to_global_vectors(family, all_dead_indices, all_percolated_indices)) {
                hout << "Error in Determine_backbone_network when calling Add_indices_to_global_vectors" << endl;
                return 0;
            }
        }
        
        //If there are GNP clusters and there are GNPs in the cluster n_cluster, find the dead GNPs
        if (HoKo->clusters_gch.size() && HoKo->clusters_gch[n_cluster].size()) {
            
            //Find the dead GNPs. Save the information on the vectors dead_gnps and percolated_gnps
            if (!Find_dead_gnps(zero_cutoff, currents_gnp, HoKo->clusters_gch[n_cluster])) {
                hout << "Error in Determine_backbone_network when calling Find_dead_gnps" << endl;
                return 0;
            }
            
            //Add the indices of percolated and non percolated CNTs
            if (!Add_gnps_to_global_vectors(family, all_dead_gnps, all_percolated_gnps)) {
                hout << "Error in Determine_backbone_network when calling Add_gnps_to_global_vectors" << endl;
                return 0;
            }
            
            /*/
            Printer *P = new Printer;
            P->Print_1d_vec(HoKo->clusters_gch[n_cluster], "cluster_gnp.txt");
            P->Print_1d_vec(percolated_gnps, "percolated_gnps.txt");
            P->Print_1d_vec(dead_gnps, "dead_gnps.txt");
            delete P;//*/
        }

    }
    
    return 1;
}
//This function determines the cutoff for "zero current"
//This step is necessary due to floating point errors
double Backbone_Network::Zero_current(const int &n_cluster, const int &R_flag, Direct_Electrifying *DEA, Hoshen_Kopelman *HoKo, const struct Electric_para &electric_param, const Cutoff_dist &cutoffs, const vector<Point_3D> &point_list, const vector<Point_3D> &point_list_gnp, const vector<double> &radii, const vector<GCH> &hybrid_particles, vector<vector<double> > &currents_cnt, vector<vector<double> > &currents_gnp)
{
    //Valiable to store the current
    double I;
    
    //Variable to store the maximum current
    double I_max = Zero;
    
    //Variable to store the cutoff for zero voltage
    double zero_cutoff;
    
    //Vector to store all the currents
    //vector<double> currents;
    
    //Check if there is any CNT cluster and the current n_cluster has CNTs
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) {
        
        //Initialize the vectors of element currents and resistances for the CNTs
        vector<double> empty_double;
        currents_cnt.assign(HoKo->clusters_cnt[n_cluster].size(), empty_double);
        
        //Calculate all currents from CNTs
        for (long int i = 0; i < (long int)HoKo->clusters_cnt[n_cluster].size(); i++) {
            
            //Current CNT
            int CNT = HoKo->clusters_cnt[n_cluster][i];
            
            //Scan all elements in the CNT
            for (long int j = 0; j < (long int)DEA->elements[CNT].size()-1; j++) {
                
                long int P1 = DEA->elements[CNT][j];
                long int P2 = DEA->elements[CNT][j+1];
                
                //Calculate the voltage difference
                //Temporarily the voltage difference will be the current value, unless the element resistance is not 1
                I = abs(Voltage_difference(P1, P2, DEA->LM_matrix, DEA->voltages));
                
                //Check the resistance flag to use the appropiate resistance
                //If R_flag is 0, the voltage difference is the current, so I remains the same
                if (R_flag) {
                    //If R_flag is 1, then use the actual resistance
                    
                    //Calculate the Resistance of the CNT segment
                    double Re = DEA->Calculate_resistance_cnt(point_list, P1, P2, radii[CNT], electric_param.resistivity_CNT);
                    //Calculate current
                    I = I/Re;
                }
                //Add current to vector of all currents
                //currents.push_back(I);
                //Check if I is the maximum current so far
                if (I > I_max) {
                    I_max = I;
                }
                //Add current to the vector of element currents
                currents_cnt[i].push_back(I);
            }
        }
    }
    
    //Check if there is any GNP cluster and the current n_cluster has GNPs
    if (HoKo->clusters_gch.size() && HoKo->clusters_gch[n_cluster].size()) {
        
        //Initialize the vectors of element currents and resistances for the GNPs
        vector<double> empty_double;
        currents_gnp.assign(HoKo->clusters_gch[n_cluster].size(), empty_double);
        
        //Calculate currents from GNPs
        for (int i = 0; i < (int)HoKo->clusters_gch[n_cluster].size(); i++) {
            
            int hyb = HoKo->clusters_gch[n_cluster][i];
            //Scan all edges of the triangulation
            for (int j = 0; j < (int)hybrid_particles[hyb].triangulation.size(); j++) {
                
                //Calculate the voltage difference
                long int P1 = hybrid_particles[hyb].triangulation[j][0];
                long int P2 = hybrid_particles[hyb].triangulation[j][1];
                
                //Determne falg of triangulation resistor
                //0: CNT-CNT tunnel
                //1: GNP-GNP tunnel
                //2: CNT-GNP tunnel
                int flag = 0, gnp_flag1 = 0, gnp_flag2 = 0;
                
                //Get the nodes
                long int node1, node2;
                //Get the points
                Point_3D point1, point2;
                //Get the radii or thickness
                double rad1, rad2;
                //flags: CNT point (1) or GNP point (0)
                if (hybrid_particles[hyb].triangulation_flags[j][0]) {
                    //Get the data from CNT
                    node1 = DEA->LM_matrix[P1];
                    point1 = point_list[P1];
                    rad1 = radii[point_list[P1].flag];
                } else {
                    //Get the data from GNP
                    node1 = DEA->LM_matrix_gnp[P1];
                    point1 = point_list_gnp[P1];
                    rad1 = hybrid_particles[hyb].gnp.hei_z;
                }
                if (hybrid_particles[hyb].triangulation_flags[j][1]) {
                    //Get the data from CNT
                    node2 = DEA->LM_matrix[P2];
                    point2 = point_list[P2];
                    rad2 = radii[point_list[P2].flag];
                } else {
                    //Get the data from GNP
                    node2 = DEA->LM_matrix_gnp[P2];
                    point2 = point_list_gnp[P2];
                    rad2 = hybrid_particles[hyb].gnp.hei_z;
                }
                
                //Calculate the voltage difference
                //Temporarily the voltage difference will be the current value, unless the element resistance is not 1
                I = abs(DEA->voltages[node2] - DEA->voltages[node1]);
                
                //Check the resistance flag to use the appropiate resistance
                //If R_flag is 0, the voltage difference is the current, so I remains the same
                if (R_flag) {
                    //If R_flag is 1, then use the actual resistance
                    
                    //Determine if the points in the GNP resistor are CNT points or GNP points
                    if (gnp_flag1 && gnp_flag2)
                        flag = 1;
                    else if (!gnp_flag1 && !gnp_flag2)
                        flag = 0;
                    else
                        flag = 2;
                    
                    //Calculate resistance of the triangulation edge
                    double Re = DEA->Calculate_resistance_gnp(flag, point1, point2, rad1, rad2, hybrid_particles[hyb], electric_param);
                    //Calculate current
                    I = I/Re;
                }
                
                //Add current to vector of all currents
                //currents.push_back(I);
                //Check if I is the maximum current so far
                if (I > I_max) {
                    I_max = I;
                }
                //Add current to the vector of element currents
                currents_gnp[i].push_back(I);
            }
        }
    }
    
    //Calculate all currents from CNT-CNT tunnels
    for (int i = 0; i < (int)DEA->elements_tunnel.size(); i++) {
        
        //Tunnel elements have only two elements per vector
        long int P1 = DEA->elements_tunnel[i][0];
        long int P2 = DEA->elements_tunnel[i][1];
        
        //Calculate the voltage difference
        //Temporarily the voltage difference will be the current value, unless the element resistance is not 1
        I = abs(Voltage_difference(P1, P2, DEA->LM_matrix, DEA->voltages));
        //Check the resistance flag to use the appropiate resistance
        //If R_flag is 0, the voltage difference is the current, so I remains the same
        if (R_flag) {
            //If R_flag is 1, then use the actual resistance
            
            //Calculate the Resistance of the tunnel
            double Re = DEA->Calculate_resistance_tunnel(radii[point_list[P1].flag], point_list[P1], radii[point_list[P2].flag], point_list[P2], electric_param, cutoffs.van_der_Waals_dist);
            
            //Calculate current
            I = I/Re;
        }
        //Add current to vector of all currents
        //currents.push_back(I);
        //Check if I is the maximum current so far
        if (I > I_max) {
            I_max = I;
        }
    }
    
    //Calculate all currents from GNP-GNP tunnels
    for (int i = 0; i < (int)DEA->elements_gnp_tunnel.size(); i++) {
        
        //Tunnel elements have only two elements per vector
        long int P1 = DEA->elements_gnp_tunnel[i][0];
        long int P2 = DEA->elements_gnp_tunnel[i][1];
        
        //Calculate the voltage difference
        //Temporarily the voltage difference will be the current value, unless the element resistance is not 1
        I = abs(Voltage_difference(P1, P2, DEA->LM_matrix_gnp, DEA->voltages));
        //Check the resistance flag to use the appropiate resistance
        //If R_flag is 0, the voltage difference is the current, so I remains the same
        if (R_flag) {
            //If R_flag is 1, then use the actual resistance
            
            //Calculate the Resistance of the tunnel
            double Re = DEA->Calculate_resistance_tunnel(1, hybrid_particles[point_list_gnp[P1].flag], point_list_gnp[P1], hybrid_particles[point_list_gnp[P2].flag].gnp.hei_z, point_list_gnp[P2], electric_param, cutoffs.van_der_Waals_dist);
            
            //Calculate current
            I = I/Re;
        }
        //Add current to vector of all currents
        //currents.push_back(I);
        //Check if I is the maximum current so far
        if (I > I_max) {
            I_max = I;
        }
    }
    
    //Calculate all currents from CNT-GNP tunnels
    for (int i = 0; i < (int)DEA->elements_mixed_tunnel.size(); i++) {
        
        //Tunnel elements have only two elements per vector
        long int P1 = DEA->elements_mixed_tunnel[i][0];
        long int P2 = DEA->elements_mixed_tunnel[i][1];
        
        //Get the nodes
        long int node1 = DEA->LM_matrix[P1];
        long int node2 = DEA->LM_matrix_gnp[P2];
        
        //Calculate the voltage difference
        //Temporarily the voltage difference will be the current value, unless the element resistance is not 1
        I = abs(DEA->voltages[node2] - DEA->voltages[node1]);
        //Check the resistance flag to use the appropiate resistance
        //If R_flag is 0, the voltage difference is the current, so I remains the same
        if (R_flag) {
            //If R_flag is 1, then use the actual resistance
            
            //Calculate the Resistance of the tunnel
            double Re = DEA->Calculate_resistance_tunnel(2, hybrid_particles[point_list_gnp[P2].flag], point_list_gnp[P2], radii[point_list[P1].flag], point_list[P1], electric_param, cutoffs.van_der_Waals_dist);
            
            //Calculate current
            I = I/Re;
        }
        //Add current to vector of all currents
        //currents.push_back(I);
        //Check if I is the maximum current so far
        if (I > I_max) {
            I_max = I;
        }
    }
    
    //Sort currents
    //sort(currents.begin(),currents.end());
    
    //The error cutoff seems to work well with a drop in 9 orders of magnitude of the current. So that is how the cutoff is set.
    //This idea comes from Li and Chou's paper of the DEA in which using a voltage of 1V, a drop in 9 orders of magnitude
    //in the current gave good results.
    //zero_cutoff = currents.back()*1e-9;
    zero_cutoff = I_max*1e-9;
    
    /*/
    Printer *P = new Printer;
    if (R_flag) {
        P->Print_1d_vec(currents, "currents_R.txt");
        P->Print_2d_vec(currents_gnp, "currents_gnp_R.txt");
    } else {
        P->Print_1d_vec(currents, "currents.txt");
        P->Print_2d_vec(currents_gnp, "currents_gnp.txt");
    }
    delete P;//*/
    
    return zero_cutoff;
}
//This function calculates the voltage difference between two nodes
double Backbone_Network::Voltage_difference(const long int &P1, const long int &P2, const vector<int> &LM_matrix, const vector<double> &voltages)
{
    //Get the node numbers
    int node1 = LM_matrix[P1];
    int node2 = LM_matrix[P2];
    //Calculate the voltage difference
    return voltages[node2] - voltages[node1];
}
//This function scans all CNTs and calculates the currents again, then it decides which CNTs are part of the backbone and which CNTs are not percolated
int Backbone_Network::Find_dead_branches(const double &zero_cutoff, const vector<vector<double> > &currents_cnt, const vector<int> &cluster, vector<vector<long int> > &elements)
{
    //This varibale is used to initialize the vectors below
    vector<long int> empty;
    //This variable will store the begining and end point of the non-conducting branches of each CNT in the cluster
    //Then each vector withing the vector can only have size 2 or 4
    dead_indices.assign(cluster.size(), empty);
    //This variable will store the begining and end point of the conducting segment of each CNT in the cluster
    //Then each vector withing the vector can only have size 2
    percolated_indices.assign(cluster.size(), empty);
    
    //Variables
    int CNT;
    //Variables used to calculate the current
    double I;
    
    //Scan all the CNTs that belong to the current cluster
    for (long int i = 0; i < (long int)cluster.size(); i++) {
        CNT = cluster[i];
        //The dead branches are in th extremes of a CNT, they cannot be in the middle.
        
        //Vector to store the elements that have a current above the cutoff
        vector<double> element_index;
        //vectors cluster and current_e have the same size
        for (int j = 0; j < (int)currents_cnt[i].size(); j++) {
            //Get the element current
            I = currents_cnt[i][j];
            if (I > zero_cutoff) {
                //If the current is above the zero cutoff, add the current j-index to the vector of element_index
                element_index.push_back(j);
            }
        }
        
        //Variable for the j-index
        int jj;
        //----------- Check for a conducting segment -----------
        if (element_index.size()) {
            //----------- Add conducting segment -----------
            //The first point of the conducting segment is given by the first j-index in element_index
            jj = element_index.front();
            percolated_indices[i].push_back(elements[CNT][jj]);
            //The second point of the conducting segment is given by the last j-index in element_index
            jj = element_index.back();
            //An element has two indices, and for the conducting segment, I need the second index of the last conducting element
            //Thus in this case I use jj+1
            percolated_indices[i].push_back(elements[CNT][jj+1]);
            
            //----------- Check for a dead branch on the front -----------
            //If the first j-index in element_index is NOT 0, then there is a dead brach at the front of the CNT
            if (element_index.front() != 0) {
                //The first index of the dead branch is the front of the CNT
                dead_indices[i].push_back(elements[CNT].front());
                //The second index is given by the first j-index in element_index
                jj = element_index.front();
                dead_indices[i].push_back(elements[CNT][jj]);
            }
            
            //----------- Check for a dead branch on the back -----------
            //If the last j-index in element_index is NOT current_e[i].size()-1, then there is a dead brach at the back of the CNT
            if (element_index.back() != (currents_cnt[i].size()-1) ) {
                //The first index of the dead branch is given by the last j-index in element_index
                jj = element_index.back();
                //The last j-index is the first index of the last conducting element
                //An element has two indices, and for the dead branch, I need the second index of the last conducting element
                //Thus in this case I use jj+1
                dead_indices[i].push_back(elements[CNT][jj+1]);
                //The last index of the dead branch is the back of the CNT
                dead_indices[i].push_back(elements[CNT].back());
            }
            
        } else {
            //----------- Whole CNT is a dead branch -----------
            //If there is no conducting segment, the whole CNT is a dead brach
            //Thus, add the first and last points of the CNT to the dead indices
            dead_indices[i].push_back(elements[CNT].front());
            dead_indices[i].push_back(elements[CNT].back());
        }
        
        //Check that the vectors have the correct size. percolated_indices can only have size 0 or 2
        if ( (percolated_indices[i].size() != 2) && (percolated_indices[i].size() != 0) ) {
            hout << "Error in Find_dead_branches. The vector percolated_indices["<<i<<"] has size " << percolated_indices[i].size();
            hout << " but it can only have size 0 or 2." << endl;
            hout << "\tThe vector dead_indices["<<i<<"] has size " << dead_indices[i].size()<<endl;
            hout << "\tCNTs in cluster: " << cluster.size() << ", current CNT: "<<CNT;
            hout << ", elements[CNT].size()=" << elements[CNT].size() << endl;
            hout << "\tLast calculated current: "<<I<<", zero cutoff: "<<zero_cutoff<<endl;
            return 0;
        }
        //Check that the vectors have the correct size. dead_indices can only have size 0, 2 or 4
        if ( (dead_indices[i].size() != 4) && (dead_indices[i].size() != 2) && (dead_indices[i].size() != 0) ) {
            hout << "Error in Find_dead_branches. The vector dead_indices["<<i<<"] has size " << dead_indices[i].size();
            hout << " but it can only have size 0, 2 or 4." <<endl;
            hout << "\tThe vector percolated_indices["<<i<<"] has size " << percolated_indices[i].size()<<endl;
            hout << "\tCNTs in cluster: " << cluster.size() << ", current CNT: "<<CNT;
            hout << ", elements[CNT].size()=" << elements[CNT].size() << endl;
            hout << "\tLast calculated current: "<<I<<", zero cutoff: "<<zero_cutoff<<endl;
            return 0;
        }
        
    }
    
    //Printer *P = new Printer;
    //P->Print_2d_vec(percolated_indices, "percolated_indices.txt");
    //P->Print_2d_vec(dead_indices, "dead_indices.txt");
    //delete P;
    
    return 1;
}
//
int Backbone_Network::Find_dead_gnps(const double &zero_cutoff, const vector<vector<double> > &currents_gnp, const vector<int> &cluster_gch)
{
    //Variables
    int hyb;
    
    //Scan every GNP in the cluster
    for (int i = 0; i < (int)cluster_gch.size() ; i++) {
        //Flag to determine if the GNP is part of the backbone or not
        int percolated = 0;
        
        //Current hybrid particle
        hyb = cluster_gch[i];
        
        //Scan every current in the resistors coming from the triangulation
        //current_gnp and cluster_gch have the same size
        for (int j = 0; j < (int)currents_gnp[i].size(); j++) {
            //If at least one current (current_gnp[i][j]) is above the cutoff then the GNP is part of the backbone
            if (currents_gnp[i][j] > zero_cutoff) {
                //Set the percolated flag to 1
                percolated = 1;
                
                //Add the GNP to the local vector of percolated GNPs
                percolated_gnps.push_back(hyb);
                //There is no need to calculate more currents so break the loop
                break;
            }
        }
        
        //Check if the percolated flag was set to 1
        if (!percolated) {
            
            //If not set to one, then add the GNP to the cluster of dead GNPs.
            dead_gnps.push_back(hyb);
        }
    }
    
    //Printer *P = new Printer;
    //P->Print_1d_vec(percolated_gnps, "percolated_gnps.txt");
    //delete P;
    
    return 1;
}
//
int Backbone_Network::Add_indices_to_global_vectors(const int &family, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices)
{
    //Add indices for percolated segments
    for (int i = 0; i < (int)percolated_indices.size(); i++) {
        for (int j = 0; j < (int)percolated_indices[i].size(); j++) {
            all_indices[family].push_back(percolated_indices[i][j]);
        }
    }
    
    //Add indices for dead branches
    for (int i = 0; i < (int)dead_indices.size(); i++) {
        for (int j = 0; j < (int)dead_indices[i].size(); j++) {
            all_dead_indices[family].push_back(dead_indices[i][j]);
        }
    }
    
    return 1;    
}
//
int Backbone_Network::Add_gnps_to_global_vectors(const int &family, vector<vector<int> > &all_dead_gnps, vector<vector<int> > &all_percolated_gnps)
{
    //Add indices for percolated GNPs
    for (int i = 0; i < (int)percolated_gnps.size(); i++) {
        all_percolated_gnps[family].push_back(percolated_gnps[i]);
    }
    
    //Add indices for dead GNPs
    for (int i = 0; i < (int)dead_gnps.size(); i++) {
        all_dead_gnps[family].push_back(dead_gnps[i]);
    }
    
    return 1;
}

