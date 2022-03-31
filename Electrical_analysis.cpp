//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Find the backbone and calculate the electrical resistivity and resistance on each direction
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Electrical_analysis.h"

int Electrical_analysis::Perform_analysis_on_clusters(const int &iter, const cuboid &window, const Simu_para &simu_param, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const Visualization_flags &vis_flags, const Output_data_flags &out_flags, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<Point_3D> &points_gnp, vector<vector<long int> > &structure_gnp, vector<GNP> &gnps)
{
    //Time variables
    time_t ct0, ct1;
    
    //Vector of parallel resistors
    //Each cluster will contribute with a resistor to each direction in which it percolates
    //So each cluster adds a parallel resistor on each percolated direction
    vector<vector<double> > parallel_resistors(3, vector<double>());
    
    //Get the number of clusters
    int n_clusters = Get_number_of_clusters(HoKo->clusters_cnt, HoKo->clusters_gnp);

    //Clear GNP triangulations (if needed)
    if (!Clear_triangulations(HoKo->clusters_gnp, gnps)) {
        hout << "Error in Perform_analysis_on_clusters when calling Clear_triangulations" << endl;
        return 0;
    }
    
    //Create a Backbone_Network object so that the vectors for nanoparticle volumes
    //and dead branches and gnps are intialized
    Backbone_Network *BN = new Backbone_Network;
    
    //Scan every percolated cluster
    for (int j = 0; j < n_clusters; j++) {
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Current iteration
        hout << "=============================" <<endl;
        hout << "\tCluster " << j+1 << " of " << n_clusters <<", family " << HoKo->family[j] << endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Direct Electrifying algorithm
        Direct_Electrifying *DEA = new Direct_Electrifying;
        
        //Resitor flag, set to 0 to use unit resistors
        int R_flag = 0;
        
        //DEA with unit resistors
        ct0 = time(NULL);
        //hout<<"DEA->Compute_voltage_field"<<endl;
        if (!DEA->Compute_voltage_field(j, R_flag, window, simu_param, electric_param, cutoffs, HoKo, Cutwins, points_cnt, radii, structure_gnp, points_gnp, gnps)) {
            hout<<"Error in Perform_analysis_on_clusters when calling DEA->Compute_voltage_field"<<endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Calculate voltage field time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //Check if triangulations are to be exported
        //Eport triangulations when there are GNPs in the current cluster and only when
        //extracting the backbone
        if (vis_flags.triangulations && gnps.size() && HoKo->clusters_gnp.size() && HoKo->clusters_gnp[j].size()) {
            
            //Export triangulations
            //hout<<"Export_triangulations"<<endl;
            if (!Export_triangulations(iter, HoKo->clusters_gnp[j], gnps, points_gnp)) {
                hout<<"Error in Compute_voltage_field when calling Export_triangulations"<<endl;
                return 0;
            }
        }
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the backbone and dead branches
        ct0 = time(NULL);
        //hout<<"BN->Determine_backbone_network"<<endl;
        if (!BN->Determine_backbone_network(j, R_flag, simu_param.simulation_scope, vis_flags.backbone, DEA->voltages, DEA->LMM_cnts, DEA->LMM_gnps, electric_param, cutoffs, structure_cnt, points_cnt, radii, points_gnp, structure_gnp, gnps, HoKo)) {
            hout<<"Error in Perform_analysis_on_clusters when calling Backbonet->Determine_backbone_network"<<endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Determine backbone network time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //Delete object to free memory
        delete DEA;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Check if it is needed to calculate the resistor network
        if (!simu_param.simulation_scope) {
            
            //Set now the R_flag to 1 to indicate that actual resistances will be used
            R_flag = 1;

            //If there is more than one cluster (i.e., j > 0), clear the triangulations
            //of the previous cluster so that it does not interfere with the current cluster
            if (j) {
                //hout << "Clear_triangulations_of_cluster" << endl;
                if (!Clear_triangulations_of_cluster(HoKo->clusters_gnp[j-1], gnps)) {
                    hout << "Error in Perform_analysis_on_clusters when calling Clear_triangulations_of_cluster" << endl;
                    return 0;
                }
            }
            
            //DEA with actual resistors along each percolated direction for current cluster
            ct0 = time(NULL);
            //hout<<"Electrical_resistance_along_each_percolated_direction"<<endl;
            if (!Electrical_resistance_along_each_percolated_direction(R_flag, j, window, HoKo, Cutwins, simu_param, electric_param, cutoffs, structure_cnt, points_cnt, radii, structure_gnp, points_gnp, gnps, parallel_resistors)) {
                hout<<"Error in Perform_analysis_on_clusters when calling Electrical_resistance_along_each_percolated_direction"<<endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Calculate resistance along each percolated direction time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        }
    }
    
    //Calculate clusters fractions
    //hout<<"Calculate_percolated_families_fractions"<<endl;
    if (!Calculate_percolated_families_fractions(out_flags.cnt_gnp_flag, structure_cnt, points_cnt, radii, gnps, HoKo, BN)) {
        hout<<"Error in Perform_analysis_on_clusters when calling Calculate_percolated_families_fractions"<<endl;
        return 0;
    }
    
    //Delete object to free memory
    delete BN;
    
    //Export visualization files for isolated particles if needed
    if (vis_flags.backbone && (HoKo->clusters_cnt.size() || HoKo->clusters_gnp.size())) {
        //hout<<"Export_isolated_particles"<<endl;
        if (!Export_isolated_particles(iter, structure_cnt, points_cnt, HoKo->isolated_cnt, gnps, HoKo->isolated_gnp)) {
            hout<<"Error in Perform_analysis_on_clusters when calling Export_isolated_particles"<<endl;
            return 0;
        }
    }
    
    //Check if it is needed to calculate the resistor network
    if (!simu_param.simulation_scope) {
        
        //Calculate the matrix resistances on each direction
        //hout<<"Calculate_resistances_and_resistivities"<<endl;
        if (!Calculate_resistances_and_resistivities(window, electric_param, parallel_resistors)) {
            hout<<"Error in Perform_analysis_on_clusters when calling Calculate_resistances_and_resistivities"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function gets the number of clusters
int Electrical_analysis::Get_number_of_clusters(const vector<vector<int> >& clusters_cnt, const vector<vector<int> >& clusters_gnp)
{
    //hout<<"clusters_cnt.size()="<<HoKo->clusters_cnt.size()<<endl;
    //hout<<"clusters_gnp()="<<HoKo->clusters_gnp.size()<<endl;
    if (clusters_cnt.size()) {
        return (int)clusters_cnt.size();

    }
    else if (clusters_gnp.size()) {
        return (int)clusters_gnp.size();
    }

    //If both vectors of clusters are empty, then there are no clusters, so return 0
    return 0;
}
//This function clears the triangulations stored in the GNP objects
int Electrical_analysis::Clear_triangulations(const vector<vector<int> >& clusters_gnp, vector<GNP>& gnps)
{
    //Clear the triangulations only if there are GNP clusters
    if (clusters_gnp.size()) 
    {
        //Iterate over all clusters
        for (size_t i = 0; i < gnps.size(); i++)
        {
            //Clear the triangulation of GNP i
            gnps[i].triangulation.clear();
        }
    }

    return 1;
}
//This function exports the tirangulation edges of a GNP cluster
int Electrical_analysis::Export_triangulations(const int &iter, const vector<int> &cluster_gnp, const vector<GNP> &gnps, const vector<Point_3D> &points_gnp)
{
    //VTK object to export visualization files
    VTK_Export VTKE;
    
    //Generate all triangulations files
    for (int i = 0; i < (int)cluster_gnp.size(); i++) {
        
        //Generate filename
        string filename = "Triangulation_" + to_string(i) + ".vtk";
        
        //Get current GNP number
        int gnp_i = cluster_gnp[i];
        
        //Export triangulation if there are edges to export
        if (gnps[gnp_i].triangulation.size()) {
            if (!VTKE.Export_triangulation(points_gnp, gnps[gnp_i].triangulation, filename)) {
                hout<<"Error in Export_tirangulations when calling VTKE.Export_triangulation, i="<<i<<endl;
                return 0;
            }
        }
    }
    
    //Make a directory for current cluster
    char command[100];
    sprintf(command, "mkdir triangulations_%.4d", iter);
    system(command);
    
    //Move all triangulation files to the corresponding folder
    sprintf(command, "mv Triangulation_*.vtk triangulations_%.4d", iter);
    system(command);
    
    return 1;
}
//This function clears the triangulations of a given cluster
int Electrical_analysis::Clear_triangulations_of_cluster(const vector<int>& cluster_gnp, vector<GNP>& gnps)
{
    //Iterate over the GNPs in the cluster
    for (size_t i = 0; i < cluster_gnp.size(); i++)
    {
        //Get GNP number
        int gnp_i = cluster_gnp[i];

        //Clear triangulation of GNP i
        gnps[gnp_i].triangulation.clear();
    }

    return 1;
}
//This function calculates the electrical resistance along each percolated direction of a cluster
int Electrical_analysis::Electrical_resistance_along_each_percolated_direction(const int &R_flag, const int &n_cluster, const cuboid &window, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const Simu_para &simu_param, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, vector<GNP> &gnps, vector<vector<double> > &paralel_resistors)
{
    //Time variables
    time_t ct0, ct1;
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Get the vector of directions
    vector<int> directions;
    if (!Vector_of_directions(HoKo->family[n_cluster], simu_param, directions)) {
        hout << "Error in Perform_analysis_on_cluster when calling Vector_of_directions" << endl;
        return 0;
    }
    
    //Save the family
    //This is actually not strictly needed, but is added in case it is needed
    //for future development that needs the original family number
    int family = HoKo->family[n_cluster];
    
    //Calculate the electrical resistance per direction
    for (int k = 0; k < (int)directions.size(); k++) {
        
        //Direct Electrifying algorithm to calculate electrical resistance
        Direct_Electrifying *DEA_Re = new Direct_Electrifying;
        
        //Temporarily set the family to be equal to percoalted direction k
        HoKo->family[n_cluster] = directions[k];
        
        //Run a new DEA to obtain the new voltage field in the backbone using the actual resistances
        //As the variable for family use the percolated direction k
        ct0 = time(NULL);
        //hout << "DEA_Re->Compute_voltage_field k=" << k << endl;
        if (!DEA_Re->Compute_voltage_field(n_cluster, R_flag, window, simu_param, electric_param, cutoffs, HoKo, Cutwins, points_cnt, radii, structure_gnp, points_gnp, gnps)) 
        {
            hout<<"Error in Electrical_resistance_along_each_percolated_direction when calling DEA_Re->Compute_voltage_field"<<endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Voltage field on backbone time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //With the new voltage field calculate the current going through a face and calculate the resistance along that direction
        //hout << "Calculate_parallel_resistor" << endl;
        if (!Calculate_parallel_resistor(directions[k], n_cluster, electric_param, DEA_Re, points_cnt, radii, HoKo->elements_cnt, HoKo->clusters_cnt, Cutwins->boundary_cnt, points_gnp, gnps, HoKo->clusters_gnp, Cutwins->boundary_gnp, paralel_resistors)) 
        {
            hout<<"Error in Perform_analysis_on_clusters when calling DEA->Calculate_parallel_resistor"<<endl;
            return 0;
        }
        //hout << "Calculate_parallel_resistor done" << endl;
        
        //Delete objects to free memory
        delete DEA_Re;
        
    }
    
    //Set the family to its original value
    //This is actually not strictly needed, but is added in case it is needed
    //for future development that needs the original family number
    HoKo->family[n_cluster] = family;
    
    return 1;
}
//This function generates a vector with the different directions in which the electrical resistance needs to be calculated
//e.g. if the family is 3 (i.e. XY), a vector with elements {0,1} is generated
int Electrical_analysis::Vector_of_directions(const int &family, const Simu_para &simu_param, vector<int> &directions)
{
    if (family == 6) {
        //If the family is 6 (XYZ); then the resistance needs to be calculated in the three directions
        //Check in which directions calculations were requested
        if (simu_param.resistances[0]) {
            directions.push_back(0);
        }
        if (simu_param.resistances[1]) {
            directions.push_back(1);
        }
        if (simu_param.resistances[2]) {
            directions.push_back(2);
        }
    } else if (family == 5) {
        //If the family is 5 (YZ); then the resistance needs to be calculated in two directions: 1 and 2
        //Check in which directions calculations were requested
        if (simu_param.resistances[1]) {
            directions.push_back(1);
        }
        if (simu_param.resistances[2]) {
            directions.push_back(2);
        }
    } else if (family == 4) {
        //If the family is 4 (XZ); then the resistance needs to be calculated in two directions: 0 and 2
        //Check in which directions calculations were requested
        if (simu_param.resistances[0]) {
            directions.push_back(0);
        }
        if (simu_param.resistances[2]) {
            directions.push_back(2);
        }
    } else if (family == 3) {
        //If the family is 5 (XY); then the resistance needs to be calculated in two directions: 0 and 1
        //Check in which directions calculations were requested
        if (simu_param.resistances[0]) {
            directions.push_back(0);
        }
        if (simu_param.resistances[1]) {
            directions.push_back(1);
        }
    } else if (family <= 2) {
        //If the family is 0, 1 or 2; then the family is the same as the only direction
        //in which resistance needs to be calculated
        //Check if resistance was requested in the direction of the cluster
        if (simu_param.resistances[family]) {
            directions.push_back(family);
        }
    } else {
        hout << "Invalid family: " << family << endl;
        return 0;
    }
    
    return 1;
}
//Calculate the resistance along a direction
//Note: direction has values 0, 1 or 2 only to represent directions X, Y and Z, respectively
//The indices of boundary_cnt are as follows:
//0,1 for X boundaries
//2,3 for Y boundaries
//4,5 for Z boundaries
//Thus 2*direction will be 0, 2 or 4, i.e. the first boundary on each direction
//Then, 2*direction+1 will be 1, 3 or 5, i.e. the second boundary on each direction
int Electrical_analysis::Calculate_parallel_resistor(const int &direction, const int &n_cluster, const Electric_para &electric_param, Direct_Electrifying *DEA, const vector<Point_3D> &points_cnt, const vector<double> &radii, vector<set<long int> > &elements, const vector<vector<int> > &clusters_cnt, const vector<vector<int> > &boundary_cnt, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const vector<vector<int> > &clusters_gnp, const vector<vector<int> > &boundary_gnp, vector<vector<double> > &paralel_resistors)
{
    //Currents from CNTs
    double I_total_cnt = 0;
    double I_total_cnt_check = 0;
    
    //Currents from GNPS
    double I_total_gnp = 0;
    double I_total_gnp_check = 0;
    
    //Get boundaries
    int b1, b2;
    if (!Get_boundaries_from_direction(direction, b1, b2)) {
        hout<<"Error in Calculate_parallel_resistor when calling Get_boundaries_from_direction"<<endl;
        return 0;
    }
    
    //Flags for outputting I_total_cnt and I_total_cnt_check
    bool flag1 = false, flag2 = false;
    
    //---------------- Currents through CNTs
    //Check if there are CNTs in the cluster
    if (clusters_cnt.size() && clusters_cnt[n_cluster].size()) {
        
        //When there are no GNPs, if any of the two boundary vectors has no CNTs there is an error as there cannot be percolation in this direction
        if (clusters_gnp.empty() && ( boundary_cnt[b1].empty() || boundary_cnt[b2].empty() ) ) {
            hout << "One boundary vector along direction " << direction << " is empty, so there cannot be percolation along that direction."<< endl;
            hout << "\t boundary_cnt[" << b1 << "].size() = " << boundary_cnt[b1].size() << endl;
            hout << "\t boundary_cnt[" << b2 << "].size() = " << boundary_cnt[b2].size() << endl;
            return 0;
        }
        
        //Calculate the current passing through boundary b1
        if (!Currents_through_boundary_cnts(electric_param, DEA, points_cnt, radii, elements, boundary_cnt[b1], I_total_cnt)) {
            hout<<"Error in Calculate_parallel_resistor when calling Currents_through_boundary_cnts (b1)"<<endl;
            return 0;
        }
        
        //=========================== CURRENT Check
        //hout <<"//=========================== CURRENT Check"<<endl;
        //Calculate the current passing through boundary b2
        if (!Currents_through_boundary_cnts(electric_param, DEA, points_cnt, radii, elements, boundary_cnt[b2], I_total_cnt_check)) {
            hout<<"Error in Calculate_parallel_resistor when calling Currents_through_boundary_cnts (b2)"<<endl;
            return 0;
        }
        
        hout << "I_cnts="<<I_total_cnt<<" direction="<<direction<<endl;
        hout << "I_cnts_check="<<I_total_cnt_check<<" direction="<<direction<<endl;
        
        //Set flag1 to true
        flag1 = true;
    }
    
    //---------------- Currents through GNPs
    //Check if there are GNPs in the cluster
    if (clusters_gnp.size() && clusters_gnp[n_cluster].size()) {
        
        //When there are no CNTs, if any of the two boundary vectors has no GNPs there is an error as there cannot be percolation in this direction
        if ( clusters_cnt.empty() && ( boundary_gnp[b1].empty() || boundary_gnp[b2].empty() )) {
            hout << "One boundary vector along direction " << direction << " is empty, so there cannot be percolation along that direction."<< endl;
            hout << "\t boundary_gnp[" << b1 << "].size() = " << boundary_gnp[b1].size() << endl;
            hout << "\t boundary_gnp[" << b2 << "].size() = " << boundary_gnp[b2].size() << endl;
            return 0;
        }
        /*hout << "boundary_gnp[b1=" << b1 << "] = ";
        for (size_t i = 0; i < boundary_gnp[b1].size(); i++)
            hout << boundary_gnp[b1][i] << " ";
        hout << endl << "boundary_gnp[b2=" << b2 << "]=";
        for (size_t i = 0; i < boundary_gnp[b2].size(); i++)
            hout << boundary_gnp[b2][i] << " ";*/
        
        //hout << "//=========================== CURRENT GNP" << endl;
        //Calculate the current passing through boundary b1, which is node 0
        if (!Currents_through_boundary_gnps(0, electric_param, DEA, points_gnp, gnps, boundary_gnp[b1], I_total_gnp)) {
            hout<<"Error in Calculate_parallel_resistor when calling Currents_through_boundary_gnps (b1)"<<endl;
            return 0;
        }
        
        //=========================== CURRENT Check GNP
        //hout <<"//=========================== CURRENT Check GNP"<<endl;
        //Calculate the current passing through boundary b2, which is node 1
        if (!Currents_through_boundary_gnps(1, electric_param, DEA, points_gnp, gnps, boundary_gnp[b2], I_total_gnp_check)) {
            hout<<"Error in Calculate_parallel_resistor when calling Currents_through_boundary_gnps (b2)"<<endl;
            return 0;
        }
        
        hout << "I_gnps="<<I_total_gnp<<" direction="<<direction<<endl;
        hout << "I_gnps_check="<<I_total_gnp_check<<" direction="<<direction<<endl;
        
        //Set flag2 to true
        flag2 = true;
    }
    

    //Calculate total currents
    double I_total = I_total_cnt + I_total_gnp;
    double I_total_check = I_total_cnt_check + I_total_gnp_check;
    
    //Only output the total currents when both CNTs and GNP currenst are calculated
    if (flag1 && flag2) {
        hout << "I_total="<<I_total<<" direction="<<direction<<endl;
        hout << "I_total_check="<<I_total_check<<" direction="<<direction<<endl;
    }
    
    //Calculate total resistance
    //To calculate the resistance in a single direction,
    //the DEA was also run using a single direction
    //thus, there are only two boundary conditions and
    //the voltage difference is equal to the input voltage
    double R_total = electric_param.applied_voltage/I_total;

    //Undo the scaling on the resistor
    //R_total = R_total / electric_param.scaling_R;
    
    //Add resistor to vector of parallel resistors
    paralel_resistors[direction].push_back(R_total);
    
    return 1;
}
//For a given direction (0 for X, 1 for Y, 2 for Z), get the two boundaries (b1 and b2)
//that need to be connected for percolation to happen along that direction
int Electrical_analysis::Get_boundaries_from_direction(const int &direction, int &b1, int &b2)
{
    //Check the direction
    if (direction == 0) {
        b1 = 2;
        b2 = 4;
    }
    else if (direction == 1) {
        b1 = 3;
        b2 = 5;
    }
    else if (direction == 2) {
        b1 = 0;
        b2 = 1;
    }
    else {
        hout<<"Error in Get_boundaries_from_direction: Invalid direction. Direction can only be 0, 1 or 2. Input was:" << direction<<endl;
        return 0;
    }
    
    return 1;
}
//This function calculates the current though CNTs at a given boundary
int Electrical_analysis::Currents_through_boundary_cnts(const Electric_para &electric_param, Direct_Electrifying *DEA, const vector<Point_3D> &points_cnt, const vector<double> &radii, vector<set<long int> > &elements, const vector<int> &boundary_cnt, double &I)
{
    //Scan all CNTs at boundary
    for (int i = 0; i < (int)boundary_cnt.size(); i++) {
        
        //Current CNT
        int CNT = boundary_cnt[i];
        //hout<<"CNT="<<CNT<<" i="<<i<<endl;
        
        //Some CNTs on the boundary might not be part of the backbone or the geometric cluster
        //First check if there are any elements on the CNT,
        //if there are no elements there, skip the CNT
        if (!elements[CNT].empty()) {
            
            //Check if the front and/or back of the CNT are in contact with the boundary
            
            //Get the points of the element at the front of the CNT
            //hout<<"elements[CNT="<<CNT<<"].size="<<elements[CNT].size()<<endl;
            set<long int>::iterator it = elements[CNT].begin();
            long int P1 = *it;
            long int P2 = *(++it);
            //hout <<"front: P1="<<P1<<" P2="<<P2<<endl;
            
            //If P1 is at a boundary with presecribed conditions, then
            //add the current of the element at the front of the CNT
            if (!Current_of_element_in_boundary(P1, P2, radii[CNT], DEA, electric_param, points_cnt, I)) {
                hout<<"Error in Currents_through_boundary_cnts when calling Current_of_element_in_boundary (front)"<<endl;
                return 0;
            }
            
            //Get the points of the element at the back of the CNT
            set<long int>::reverse_iterator rit = elements[CNT].rbegin();
            P1 = *rit;
            P2 = *(++rit);
            //hout <<"back: P1="<<P1<<" P2="<<P2<<endl;
            
            //Add the current of the element at the back of the CNT
            if (!Current_of_element_in_boundary(P1, P2, radii[CNT], DEA, electric_param, points_cnt, I)) {
                hout<<"Error in Currents_through_boundary_cnts when calling Current_of_element_in_boundary (front)"<<endl;
                return 0;
            }
        }
    }
    
    return 1;
}
//This function checks if an element is at a boundary, and if so it calculates the current
//It is assumed that P1 is the point that is at either the front or the back of the CNT, only these points can be in contact with the boundary
int Electrical_analysis::Current_of_element_in_boundary(const long int &P1, const long int &P2, const double &radius, Direct_Electrifying *DEA, const Electric_para &electric_param, const vector<Point_3D> &points_cnt, double &I)
{
    //Get the node number of the first point, if the point is in the LMM matrix
    map<long int, long int>::const_iterator it = DEA->LMM_cnts.find(P1);
    if (it == DEA->LMM_cnts.end()) {
        //P1 is not in the LMM matrix, then terminate the function
        return 1;
    }
    long int node1 = DEA->LMM_cnts.at(P1);
    //hout<<"node1="<<node1<<endl;
    
    //Check where is P1
    //If P1 is node 0 or 1, then the current is calculated
    //This means it is on a valid boundary
    if (node1 <= 1) {
        
        //hout<<"P1="<<P1<<", "<<points_cnt[P1].str()<<" CNT1="<<points_cnt[P1].flag<<" node1="<<node1<<endl;
        //Get the node number of the second point, if the point is in the LMM matrix
        it = DEA->LMM_cnts.find(P2);
        if (it == DEA->LMM_cnts.end()) {
            //hout<<endl;
            //P2 is not in the LMM matrix, then terminate the function
            return 1;
        }
        long int node2 = DEA->LMM_cnts.at(P2);
        //hout<<" node2="<<node2<<endl;
        //hout<<" P2="<<P2<<", "<<points_cnt[P2].str()<<" CNT2="<<points_cnt[P2].flag<<" node2="<<node2<<endl;
        
        //Calculate voltage difference on the element,
        //In the first calculation of current, node1 is at the boundary with voltage 0
        //so the voltage drop is from node2 to node1
        double V = DEA->voltages[node2] - DEA->voltages[node1];
        
        //hout << "V1="<<DEA->voltages[node1]<<" V2="<<DEA->voltages[node2]<<"\nDV=" << V<<endl;
        //Calculate resistance of the element
        double Re;
        if (P1 > P2) {
            //Check if the resistor is at the back of the CNT, since in that case the calculated resistance will be zero;
            //this because if P1 > P2, the function that calculates the resistance does not find points after P1
            //When the resistor is at the back of the element P1 > P2, so in that case invert the point numbers
            if (!DEA->Calculate_resistance_cnt(1, points_cnt, P2, P1, radius, electric_param, Re)) {
                hout<<"Error in Current_of_element_in_boundary when calling DEA->Calculate_resistance_cnt (1)"<<endl;
                return 0;
            }
        }
        else {
            if (!DEA->Calculate_resistance_cnt(1, points_cnt, P1, P2, radius, electric_param, Re)) {
                hout<<"Error in Current_of_element_in_boundary when calling DEA->Calculate_resistance_cnt (2)"<<endl;
                return 0;
            }
        }
        
        //Calculate current and add it to the total current
        //hout<<"Re="<<Re<<" I=("<<DEA->voltages[node2]<<"-"<<DEA->voltages[node1]<<")/Re=\t"<<V/Re<<endl;
        I = I + V/Re;
    }
    
    return 1;
}
//This function calculates the current though GNPs at a given boundary
int Electrical_analysis::Currents_through_boundary_gnps(const long int &node, const Electric_para &electric_param, Direct_Electrifying *DEA, const vector<Point_3D> &points_gnp, const vector<GNP> &gnps, const vector<int> &boundary_gnp, double &I)
{
    //Iterate over the GNPs at the boundary
    for (int i = 0; i < (int)boundary_gnp.size(); i++) {
        
        //Current GNP
        int GNPi = boundary_gnp[i];
        //hout <<endl<< "GNPi=" << GNPi << " triangulation.size=" << gnps[GNPi].triangulation.size() << endl;
        
        //Some GNPs on the boundary might not be part of the backbone or the geometric cluster
        //First check if there are any triangulation edges on the GNP,
        //if there are no edges there, skip the GNP
        if (gnps[GNPi].triangulation.size()) {
            
            //Iterate over the triangulation edges
            for (int j = 0; j < (int)gnps[GNPi].triangulation.size(); j++) {
                
                //Get the vertices of the triangulation
                long int v1 = gnps[GNPi].triangulation[j].v1;
                long int v2 = gnps[GNPi].triangulation[j].v2;
                //hout << "v1=" << v1 << " v2=" << v2 << endl;
                
                //Get the nodes of the vertices of the triangulation
                long int nodeA = DEA->LMM_gnps.at(v1);
                //hout<<"nodeA="<<nodeA<<" v1="<<v1<<" P(v1)="<<points_gnp[v1].str()<<endl;
                long int nodeB = DEA->LMM_gnps.at(v2);
                //hout<<"nodeB="<<nodeB<<" v2="<<v2<<" P(v2)="<<points_gnp[v2].str()<<endl;
                
                //Check if any node is at a boundary
                if (nodeA <= 1 || nodeB <= 1) {
                    
                    //Check that the nodes are not in the same boundary
                    if (nodeA == nodeB) {
                        hout<<"Error in Currents_through_boundary_gnps, nodes are at the same boundary: node1="<<nodeA<<", node2="<<nodeB<<endl;
                        return 0;
                    }
                    
                    //Set node1 to be the node at boundary b
                    long int node1 = (nodeA == node)? nodeA: nodeB;
                    long int node2 = (nodeA == node)? nodeB: nodeA;
                    //hout<<"node1="<<node1<<" node2="<<node2<<endl;
                    
                    //Get the radii for calculating the resistance of the triangulation edge
                    double rad1 = (DEA->points_cnt_rad.find(v1) == DEA->points_cnt_rad.end())? gnps[GNPi].t/2: DEA->points_cnt_rad.at(v1);
                    double rad2 = (DEA->points_cnt_rad.find(v2) == DEA->points_cnt_rad.end())? gnps[GNPi].t/2: DEA->points_cnt_rad.at(v2);
                    
                    //Calculate the drop in voltage
                    //In the first calculation of current, node1 is at the boundary with voltage 0
                    //so the voltage drop is from node2 to node1
                    double V = DEA->voltages[node2] - DEA->voltages[node1];
                    
                    //Calculate resistance of triangulation edge
                    double Re;
                    if (!DEA->Calculate_resistance_gnp(points_gnp[v1], points_gnp[v2], rad1, rad2, electric_param, Re)) {
                        hout<<"Error in Currents_through_boundary_gnps when calling DEA->Calculate_resistance_gnp"<<endl;
                        return 0;
                    }
                    //hout<<"V="<<DEA->voltages[node2]<<" - "<<DEA->voltages[node1]<<", Re="<<Re<<" L="<< points_gnp[v1].distance_to(points_gnp[v2])<< endl;
                    
                    //Add calcualted current
                    I = I + V/Re;
                    //hout<<"I="<<I<<endl;
                }
            }
        }
    }
    
    return 1;
}
//This function calculates the fractions for each percoalted family,
//which is used to calculat the effiency
int Electrical_analysis::Calculate_percolated_families_fractions(const int &cnt_gnp_flag, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<GNP> &gnps, Hoshen_Kopelman *HoKo, Backbone_Network *BN)
{
    //Total volume of dead CNTs
    double total_dead_branches = 0;
    
    //Total volume od dead GNPs
    double total_dead_gnps = 0;
    
    //Add the volumes of dead branches and GNPs from percoalted clusters
    for (int i = 0; i < (int)BN->dead_branches.size(); i++) {
        total_dead_branches = total_dead_branches + BN->dead_branches[i];
        total_dead_gnps = total_dead_gnps + BN->dead_gnps[i];
    }
    
    //Calculate the volumes of isolated CNTs and non-percoalted clusters
    double np_cnts = 0;
    if (!Calculate_volume_of_non_percolated_cnts(structure_cnt, points_cnt, radii, HoKo->isolated_cnt, BN, np_cnts)) {
        hout<<"Error in Calculate_percolated_families_fractions while calling Calculate_volume_of_non_percolated_cnts"<<endl;
        return 0;
    }
    
    //Add all the volumes of dead particles to the last element of CNT volumes
    BN->volumes_cnt.back() = total_dead_branches + np_cnts;
    
    //Calculate the volumes of isolated GNPs and non-percoalted clusters
    double np_gnps = 0;
    if (!Calculate_volume_of_non_percolated_gnps(gnps, HoKo->isolated_gnp, np_gnps)) {
        hout<<"Error in Calculate_percolated_families_fractions while calling Calculate_volume_of_non_percolated_gnps"<<endl;
        return 0;
    }
    
    //Add all the volumes of dead particles to the last element of CNT volumes
    BN->volumes_gnp.back() = total_dead_gnps + np_gnps;
    
    //Vector with volumes of percoalted particles
    vector<double> total_volumes(8,0);
    
    //Total volume of particles
    double total_volume = 0;
    double total_cnts = 0, total_gnps = 0;
    
    //Add the volumes of percolated particles
    for (int i = 0; i < (int)BN->volumes_cnt.size(); i++) {
        total_volumes[i] = total_volumes[i] + BN->volumes_cnt[i] + BN->volumes_gnp[i];
        total_volume = total_volume + BN->volumes_cnt[i] + BN->volumes_gnp[i];
        total_cnts = total_cnts + BN->volumes_cnt[i];
        total_gnps = total_gnps + BN->volumes_gnp[i];
    }
    
    //Calculate fractions
    vector<double> fractions(8,0);
    vector<double> fractions_cnts(8,0), fractions_gnps(8,0);
    for (int i = 0; i < (int)fractions.size(); i++) {
        fractions[i] = total_volumes[i]/total_volume;
        fractions_cnts[i] = BN->volumes_cnt[i]/total_cnts;
        fractions_gnps[i] = BN->volumes_gnp[i]/total_gnps;
    }
    
    //Printer object to write to files
    Printer P;
    
    //Output the volumes of percolated families if they are requested (through cnt_gnp_flag)
    if (cnt_gnp_flag && HoKo->clusters_cnt.size() && HoKo->clusters_gnp.size()) {
        
        //Print CNT and GNP volumes
        P.Append_1d_vec(BN->volumes_cnt, "volumes_cnt.txt");
        P.Append_1d_vec(BN->volumes_gnp, "volumes_gnp.txt");
        
        //Print CNT and GNP fractions
        P.Append_1d_vec(fractions_cnts, "fractions_cnt.txt");
        P.Append_1d_vec(fractions_gnps, "fractions_gnp.txt");
    }
    
    //Print volumes and fractions
    P.Append_1d_vec(total_volumes, "volumes.txt");
    P.Append_1d_vec(fractions, "fractions.txt");
    
    return 1;
}
//This function calcualtes the volumes of the CNTs that are isolated or in non-percolated clusters
int Electrical_analysis::Calculate_volume_of_non_percolated_cnts(const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<int> > &isolated_cnts, Backbone_Network *BN, double &np_cnts)
{
    //Iterate over the clusters of non-percolated CNTs
    for (int i = 0; i < (int)isolated_cnts.size(); i++) {
        
        //Iterate over the CNTs in cluster i
        for (int j = 0; j < (int)isolated_cnts[i].size(); j++) {
            
            //Get the current CNT
            int CNTj = isolated_cnts[i][j];
            
            //Get the two endpoints of the CNT
            long int P1 = structure_cnt[CNTj].front();
            long int P2 = structure_cnt[CNTj].back();
            
            //Calculate the volume of the CNT
            double cnt_vol = 0;
            if (!BN->CNT_volume_between_two_points(P1, P2, radii[CNTj], points_cnt, cnt_vol)) {
                hout<<"Error in Calculate_volume_of_non_percolated_cnts while calling BN->CNT_volume_between_two_points"<<endl;
                return 0;
            }
            
            //Add the calculated volume to the total volume
            np_cnts = np_cnts + cnt_vol;
        }
    }
    
    return 1;
}
//This function calcualtes the volumes of the GNPs that are isolated or in non-percolated clusters
int Electrical_analysis::Calculate_volume_of_non_percolated_gnps(const vector<GNP> &gnps, const vector<vector<int> > &isolated_gnps, double &np_gnps)
{
    //Iterate over the clusters of non-percolated CNTs
    for (int i = 0; i < (int)isolated_gnps.size(); i++) {
        
        //Iterate over the GNPs in cluster i
        for (int j = 0; j < (int)isolated_gnps[i].size(); j++) {
            
            //Get current GNP
            int GNPj = isolated_gnps[i][j];
            
            //Add the volume of the GNP to the total volume
            np_gnps = np_gnps + gnps[GNPj].volume;
        }
    }
    
    return 1;
}
//This function prepares the vectors needed to export isolated particles
int Electrical_analysis::Export_isolated_particles(const int &iter, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<vector<int> > &isolated_cnts, const vector<GNP> &gnps, const vector<vector<int> > &isolated_gnps)
{
    //VTK object to export visualization files
    VTK_Export VTK_E;
    
    //Prepare filenames
    string str_cnts = "isolated_cnts_from_backbone.vtk";
    string str_gnps = "isolated_gnps_from_backbone.vtk";
    
    //Create a vector with indices of all isolated CNTs
    vector<vector<long int> > all_isolated_cnts;
    //hout<<"all_isolated_cnts"<<endl;
    for (int i = 0; i < (int)isolated_cnts.size(); i++) {
        for (int j = 0; j < (int)isolated_cnts[i].size(); j++) {
            
            //Get CNT
            int CNTij = isolated_cnts[i][j];
            
            //Create a vector with the two endpoints of the CNT
            vector<long int> tmp(2);
            tmp[0] = structure_cnt[CNTij].front();
            tmp[1] = structure_cnt[CNTij].back();
            all_isolated_cnts.push_back(tmp);
        }
    }
    
    //Export the isolated CNTs using the indices (if there are CNTs)
    if (isolated_cnts.size())
    {
        //hout<<"VTK_E.Export_from_cnt_indices"<<endl;
        if (!VTK_E.Export_from_cnt_indices(points_cnt, all_isolated_cnts, str_cnts)) {
            hout << "Error in Export_isolated_particles when calling VTK_E.Export_from_cnt_indices" << endl;
            return 0;
        }
    }
    
    //Create a vector with indices of all isolated GNPs
    vector<int> all_isolated_gnps;
    //hout<<"all_isolated_gnps"<<endl;
    for (int i = 0; i < (int)isolated_gnps.size(); i++) {
        for (int j = 0; j < (int)isolated_gnps[i].size(); j++) {
            all_isolated_gnps.push_back(isolated_gnps[i][j]);
        }
    }
    
    //Export the isolated GNPs as a cluster (if there are GNPs)
    if (isolated_gnps.size())
    {
        //hout<<"VTK_E.Export_gnps_in_cluster"<<endl;
        if (!VTK_E.Export_gnps_in_cluster(gnps, all_isolated_gnps, str_gnps)) {
            hout << "Error in Export_isolated_particles when calling VTK_E.Export_gnps_in_cluster" << endl;
            return 0;
        }
    }
    
    //Move all visualization files to the folder of the iteration
    //hout<<"system commands"<<endl;
    char command[100];
    sprintf(command, "mkdir backbone_%.4d", iter);
    system(command);
    sprintf(command, "mv backbone*.vtk dead*.vtk isolated*.vtk backbone_%.4d", iter);
    system(command);
    
    return 1;
}
//This function calculates the resistance on each direction from the vector of parallel resistors
int Electrical_analysis::Calculate_resistances_and_resistivities(const cuboid &window, const Electric_para &electric_param, const vector<vector<double> > &paralel_resistors)
{
    //Variables to store the calculcated resistances and resistivities of the sample along
    //each direction
    vector<double> resistors(3,0);
    vector<double> resistivities(3,0);
    
    //Scan each direction
    //Note that matrix_resistances and paralel_resistors have the same size
    for (int i = 0; i < (int)paralel_resistors.size(); i++) {
        if (paralel_resistors[i].empty()) {
            
            //There are no resistors in direction i
            //The resistivity is the same as the matrix
            resistivities[i] = electric_param.resistivity_matrix;
            
            //Calculate resistance of the matrix along direction i
            if (!Calculate_matrix_resistance(i, window, electric_param.resistivity_matrix, resistors[i])) {
                hout<<"Error in Calculate_resistances_and_resistivities when calling Calculate_matrix_resistance"<<endl;
                return 0;
            }
        } else if (paralel_resistors[i].size() == 1) {
            
            //If direction i has only one resistor, then that is the resistance in that direction
            resistors[i] = paralel_resistors[i].front();
            
            //Calculate the resistivity
            if (!Calculate_resistivity(i, window, resistors[i], resistivities[i])) {
                hout<<"Error in Calculate_resistances_and_resistivities when calling Calculate_resistivity (1 resistor)"<<endl;
                return 0;
            }
        } else {
            //If there is more than one resistor, then calculate the
            //equivalent resistance of parallel resistors
            double R = 0;
            for (int j = 0; j < (int)paralel_resistors[i].size(); j++) {
                R = R + 1/paralel_resistors[i][j];
            }
            resistors[i] = 1/R;
            
            //Calculate the resistivity
            if (!Calculate_resistivity(i, window, resistors[i], resistivities[i])) {
                hout<<"Error in Calculate_resistances_and_resistivities when calling Calculate_resistivity (1 resistor)"<<endl;
                return 0;
            }
        }
    }
    
    //Print resistances and resistivities into a file
    Printer P;
    P.Append_1d_vec(resistors, "resistances.txt");
    P.Append_1d_vec(resistivities, "resistivities.txt");
    
    return 1;
}
//This function calculates the matrix resistance on a given direction
int Electrical_analysis::Calculate_matrix_resistance(const int &direction, const cuboid &window, const double &matrix_resistivity, double &R_M)
{
    //Lambda function to calculate the resistance of the matrix in a given direction
    auto res_matrix = [](const double &length, const double &A, const double &rho_m){
        return (rho_m*length/A);
    };
    
    //------------------ Resistance along the x-direction
    if (direction == 0) {
        R_M = res_matrix(window.len_x, window.wid_y*window.hei_z, matrix_resistivity);
    }
    //------------------ Resistance along the y-direction
    else if (direction == 1) {
        R_M = res_matrix(window.wid_y, window.len_x*window.hei_z, matrix_resistivity);
    }
    //------------------ Resistance along the z-direction
    else if (direction == 2) {
        R_M = res_matrix(window.hei_z, window.len_x*window.wid_y, matrix_resistivity);
    }
    else {
        hout<<"Error in Calculate_matrix_resistance: invalid direction. Direction can only be 0, 1, or 2. Input was"<<direction<<endl;
        return 0;
    }
    
    return 1;
}
//This function calculates resistivity on a given direction
int Electrical_analysis::Calculate_resistivity(const int &direction, const cuboid &window, const double &resistance, double &rho)
{
    //Lambda function to calculate the resistivity in a given direction
    auto rho_dir = [](const double &length, const double &A, const double &R){
        return (R*A/length);
    };
    
    //------------------ Resistivity along the x-direction
    if (direction == 0) {
        rho = rho_dir(window.len_x, window.wid_y*window.hei_z, resistance);
    }
    //------------------ Resistivity along the y-direction
    else if (direction == 1) {
        rho = rho_dir(window.wid_y, window.len_x*window.hei_z, resistance);
    }
    //------------------ Resistivity along the z-direction
    else if (direction == 2) {
        rho = rho_dir(window.hei_z, window.len_x*window.wid_y, resistance);
    }
    else {
        hout<<"Error in Calculate_resistivity: invalid direction. Direction can only be 0, 1, or 2. Input was"<<direction<<endl;
        return 0;
    }
    
    
    return 1;
}
