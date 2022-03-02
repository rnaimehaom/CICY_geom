//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Implementation of Hoshen-Kopelman Algorithm
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Hoshen_Kopelman.h"

/*
 
 This function implements the Hoshen-Kopelman algorithm and generates the vector of points that are in contact.
 It also creates a connectivity vector of the points in contact. This connectivity vector is used to define the elements used in the direct electrifying algorithm.
 */

//This function groups nanoparticles into clusters
int Hoshen_Kopelman::Determine_clusters_and_percolation(const int &iter, const cuboid& sample, const Simu_para &simu_para, const Cutoff_dist &cutoffs, const Visualization_flags &vis_flags, const vector<int> &cnts_inside, const vector<vector<long int> > &sectioned_domain_cnt, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<int> > &boundary_cnt, const vector<int> &gnps_inside, const vector<vector<int> > &sectioned_domain_gnp, const vector<GNP> &gnps, const vector<vector<int> > &boundary_gnp, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp)
{
    //Label vectors for CNTs and GNPs
    vector<int> labels_cnt, labels_gnp;
    
    //Variable to count the labels from CNT clusters
    int n_labels_cnt = 0;
    
    //Label to store the total number of labels in case there are both CNTs and GNPs
    int n_total_labels = 0;
    
    //Adjacency matrix of labels for merging CNT and GNP labels
    vector<set<int> > adj_labels;
    
    //Variable to store the total number of clusters
    int n_clusters = 0;
    
    //Check if there are CNTs in the sample and make the CNT clusters accordingly
    if (simu_para.particle_type != "GNP_cuboids") {
        
        /* /If hybrid particles are added, a pre-processing is needed before making the CNT clusters
        if (simu_para.particle_type == "Hybrid_particles") {
            //
        }*/
        
        //The size of the vector labels has to be equal to the number of CNTs
        //-1 is the unassigned value for a label
        labels_cnt.assign(structure_cnt.size(), -1);
        
        //There are CNTs, make clusters
        //hout<<"Make_cnt_clusters"<<endl;
        if (!Make_cnt_clusters(points_cnt, radii, cutoffs, sectioned_domain_cnt, structure_cnt, labels_cnt, n_labels_cnt)) {
            hout<<"Error in Determine_clusters when calling Make_cnt_clusters"<<endl;
            return 0;
        }
        
        //Provisionally, the total number of labels is the same as the number of CNT labels
        //This will change if there are GNPs
        n_total_labels = n_labels_cnt;
        
        //Also, provisonally the number of clusters is the same as the number of CNT labels
        //If GNPs are generated, n_clusters will be updated in Make_mixed_clusters
        n_clusters = n_labels_cnt;
    }
    
    //Check if there are GNPs in the sample and make the GNP clusters accordingly
    if (simu_para.particle_type != "CNT_wires" && simu_para.particle_type != "CNT_deposit") {
        
        //Vectors of labels
        //The size of the vector labels has to be equal to the number of GNPs
        //-1 is the unassigned value for a label
        labels_gnp.assign(gnps.size(), -1);
        
        //Make GNP clusters
        //hout<<"Make_gnp_clusters n_total_labels="<<n_total_labels<<endl;
        if (!Make_gnp_clusters(sample, gnps_inside, sectioned_domain_gnp, gnps, cutoffs, n_labels_cnt, n_total_labels, structure_gnp, points_gnp, labels_gnp)) {
            hout<<"Error in Determine_clusters when calling Make_gnp_clusters"<<endl;
            return 0;
        }
        
        //If there are no CNTs, then the number of clusters is the same as n_total_labels
        //If CNTs were generated, n_clusters will be updated in Make_mixed_clusters
        //hout<<"clusters_gnp.size()="<<clusters_gnp.size()<<" n_total_labels="<<n_total_labels<<endl;
        n_clusters = n_total_labels;
    }
    
    //Check if there are both CNTs and GNPs in the sample and merge CNT and GNP clusters as needed
    if (simu_para.particle_type == "GNP_CNT_mix" || simu_para.particle_type == "Hybrid_particles") {
        
        //Make mixed clusters
        //hout<<"Make_mixed_clusters"<<endl;
        if (!Make_mixed_clusters(n_total_labels, cutoffs, points_cnt, radii, sectioned_domain_cnt, gnps, sectioned_domain_gnp, labels_cnt, labels_gnp, structure_gnp, points_gnp, n_clusters)) {
            hout<<"Error in Determine_clusters when calling Make_mixed_clusters"<<endl;
            return 0;
        }
    }
    
    //Variables to store the clusters
    vector<vector<int> > clusters_cnt_tmp, clusters_gnp_tmp;
    
    //Generate cluster vectors
    if (labels_cnt.size()) {
        //hout<<"Make_particle_clusters CNTs"<<endl;
        if (!Make_particle_clusters(n_clusters, cnts_inside, labels_cnt, isolated_cnt, clusters_cnt_tmp)) {
            hout<<"Error in Determine_clusters when calling Make_particle_clusters (CNT)"<<endl;
            return 0;
        }
    }
    if (labels_gnp.size()) {
        //hout<<"Make_particle_clusters GNPs"<<endl;
        if (!Make_particle_clusters(n_clusters, gnps_inside, labels_gnp, isolated_gnp, clusters_gnp_tmp)) {
            hout<<"Error in Determine_clusters when calling Make_particle_clusters (GNP)"<<endl;
            return 0;
        }
    }
    
    //Export visualization files for clusters if needed
    if (vis_flags.clusters) {
        if (!Export_clusters(0, iter, clusters_cnt_tmp, structure_cnt, points_cnt, clusters_gnp_tmp, gnps)) {
            hout<<"Error in Determine_clusters when calling Export_clusters (before percolation)"<<endl;
            return 0;
        }
    }
    
    //Determine percolation
    //hout<<"Find_percolated_clusters"<<endl;
    map<int,int> percolated_labels;
    if (!Find_percolated_clusters(n_clusters, boundary_cnt, boundary_gnp, labels_cnt, labels_gnp, clusters_cnt_tmp, clusters_gnp_tmp, percolated_labels)) {
        hout<<"Error in Determine_clusters when calling Find_percolated_clusters"<<endl;
        return 0;
    }
    
    //Export visualization files for percolated clusters if needed
    if (vis_flags.percolated_clusters) {
        if (!Export_clusters(1, iter, clusters_cnt, structure_cnt, points_cnt, clusters_gnp, gnps)) {
            hout<<"Error in Determine_clusters when calling Export_clusters (after percolation)"<<endl;
            return 0;
        }
    }
    
    //Group junctions if there is percolation
    if ( clusters_cnt.size() + clusters_gnp.size() ) {
        
        //There is at least one percoalted cluster
        //Group the junctions into the clusters they belong to
        //hout<<"Group_junctions"<<endl;
        hout << "There are percolated clusters" << endl;
        if (!Group_junctions(points_cnt, points_gnp, labels_cnt, labels_gnp, percolated_labels)) {
            hout<<"Error in Find_percolated_clusters when calling Group_junctions"<<endl;
            return 0;
        }
    }
    else {
        
        //There is no percolation, so output a message
        hout<<"There are no percolated clusters"<<endl;
    }
    
    return 1;
}
//This functions makes CNT clusters
int Hoshen_Kopelman::Make_cnt_clusters(const vector<Point_3D> &points_cnt, const vector<double> &radii, const Cutoff_dist &cutoffs, const vector<vector<long int> > &sectioned_domain_cnt, const vector<vector<long int> > &structure_cnt, vector<int> &labels_cnt, int &n_labels_cnt)
{
    //Vector of labels of labels
    vector<int> labels_labels_cnt;
    
    //Variables to compress contact segments
    vector<map<int, set<long int> > > contact_elements(radii.size());
    map<long int, map<int, long int> > point_contacts;
    map<long int, map<int, double> > point_contacts_dist;
    
    //Label the CNTs
    //hout<<"Label_cnts_in_window"<<endl;
    if (!Label_cnts_in_window(points_cnt, radii, cutoffs, sectioned_domain_cnt, point_contacts, point_contacts_dist, contact_elements, labels_cnt, labels_labels_cnt)) {
        hout<<"Error in Make_cnt_clusters when calling Label_cnts_in_window"<<endl;
        return 0;
    }
    
    //Clean up the labels to find the proper labels, i.e. merged and consecutive labels starting at 0
    //hout<<"Cleanup_labels"<<endl;
    if (!Cleanup_labels(labels_labels_cnt, labels_cnt, n_labels_cnt)) {
        hout << "Error in Make_cnt_clusters when calling Cleanup_labels" << endl;
        return 0;
    }
    
    //Compress CNT-CNT contact segments and add junctions
    //hout<<"Compress_cnt_cnt_contact_segments"<<endl;
    if (!Compress_cnt_cnt_contact_segments(cutoffs, point_contacts, point_contacts_dist, contact_elements)) {
        hout << "Error in Make_cnt_clusters when calling Compress_cnt_cnt_contact_segments" << endl;
        return 0;
    }
    
    //Add first and last points to the elements
    //hout<<"Complete_cnt_elements"<<endl;
    if (!Complete_cnt_elements(structure_cnt)) {
        hout << "Error in Make_cnt_clusters when calling Complete_cnt_elements" << endl;
        return 0;
    }
    
    return 1;
}
//This function labels the CNTs in a window
//These labels are used to make clusters
int Hoshen_Kopelman::Label_cnts_in_window(const vector<Point_3D> &points_cnt, const vector<double> &radii, const Cutoff_dist& cutoffs, const vector<vector<long int> > &sectioned_domain_cnt, map<long int, map<int, long int> > &point_contacts, map<long int, map<int, double> > &point_contacts_dist, vector<map<int, set<long int> > > &contact_elements, vector<int> &labels_cnt, vector<int> &labels_labels_cnt)
{
    //new_label will take the value of the newest cluster
    int new_label = 0;
    
    //Scan every overlapping sub-region
    for (int i = 0; i < (int)sectioned_domain_cnt.size(); i++) {
        
        //Size of subregion i
        int inner = (int)sectioned_domain_cnt[i].size();
        
        //Scan all points in subregion i
        for (int j = 0; j < inner-1; j++) {
            
            //Get first point and its CNT number
            long int P1 = sectioned_domain_cnt[i][j];
            int CNT1 = points_cnt[P1].flag;
            //hout<<"P1="<<points_cnt[P1].str()<<" sectioned_domain_cnt[i][j]="<<sectioned_domain_cnt[i][j]<<endl;
            
            //For each point P1=sectioned_domain[i][j] scan all remaning points in sectioned_domain[i]
            for (int k = j+1; k<inner; k++) {
                
                //Get second point and its CNT number
                long int P2 = sectioned_domain_cnt[i][k];
                int CNT2 = points_cnt[P2].flag;
                //hout<<"P2="<<points_cnt[P2].str()<<" sectioned_domain_cnt[i][k]="<<sectioned_domain_cnt[i][k]<<endl;
                //hout <<"P1="<<P1<<" CNT1="<<CNT1<<" P2="<<P2<<" CNT2="<<CNT2;
                //hout<<" r1="<<radii[CNT1]<<" r2="<<radii[CNT2]<<endl;
                
                //Check if:
                //The points belong to different CNTs, only when the CNTs are different
                //the distance between points is calculated
                //AND
                //The contact CNT1-CNT2 has not been visited yet
                if (CNT1 != CNT2) {
                    
                    //Calculate the distance between points
                    double dist_cnts = points_cnt[P1].distance_to(points_cnt[P2]);
                    //hout <<" dist_junc="<<dist_junc<<endl;

                    //Calculate the junction distance
                    double dist_junc = dist_cnts - radii[CNT1] - radii[CNT2];
                    if (dist_junc < cutoffs.van_der_Waals_dist)
                    {/* /
                        hout << "Junction distance too small dist_junc=" << dist_junc << " dist_cnts=" << dist_cnts << endl;
                        hout<<"P1="<<points_cnt[P1].str()<<" sectioned_domain_cnt[i="<<i<<"][j="<<j<<"]="<<sectioned_domain_cnt[i][j]<<endl;
                        hout<<"P2="<<points_cnt[P2].str()<<" sectioned_domain_cnt[i="<<i<<"][k="<<k<<"]="<<sectioned_domain_cnt[i][k]<<endl;
                        hout <<"P1="<<P1<<" CNT1="<<CNT1<<" P2="<<P2<<" CNT2="<<CNT2;
                        hout<<" r1="<<radii[CNT1]<<" r2="<<radii[CNT2]<<endl;
                        return 0;// */

                        //The distance between points might be less than the van der Waals
                        //distance if penetrations are allowed or Abaqus causes them
                        dist_junc = cutoffs.van_der_Waals_dist;
                    }
                    
                    //Compare the distance between points against the junction distance
                    //to check whether there is tunneling
                    if (dist_junc <= cutoffs.tunneling_dist) {
                        
                        //Here is where the actual HK76 algorithm takes place
                        if (!HK76(CNT1, CNT2, new_label, labels_cnt, labels_labels_cnt)) {
                            hout << "Error in Label_cnts_in_window when calling HK76" << endl;
                            return 0;
                        }
                        
                        //Sort the CNTs as maximum and minimum
                        int CNTa = min(CNT1, CNT2);
                        long int Pa = (CNT1 < CNT2)? P1: P2;
                        int CNTb = max(CNT1, CNT2);
                        long int Pb = (CNT2 < CNT1)? P1: P2;
                        //hout <<"CNT1="<<CNT1<<" P1="<<P1<<" CNT2="<<CNT2<<" P2="<<P2<<endl;
                        //hout<<"CNTa="<<CNTa<<" Pa="<<Pa<<" CNTb="<<CNTb<<" Pb="<<Pb<<endl;
                        
                        //Add points to the elements used for contacts
                        //Only add the points to the CNT with lowest number
                        contact_elements[CNTa][CNTb].insert(Pa);
                        
                        
                        //Make sure whether the contact has been added to the contact_points map
                        //and keep the smallest separation
                        //This map eliminates the posibility of a point having two contacts
                        //with the same CNT
                        //Also, this map is used to compress "contact segments"
                        //
                        //Clauses:
                        //A: point_contacts.find(Pa) == point_contacts.end()
                        //B: point_contacts[Pa].find(CNTb) == point_contacts[P1].end()
                        //C: point_contacts_dist[Pa][CNTb] < junction_dist
                        //
                        //If A is true, then Pa is not in point_contacts, thus the statements inside the
                        //if-statement are evaluated and Pa is added to point_contacts and CNTb is added
                        //to point_contacts_dist[Pa]
                        //
                        //If A, is not true, then Pa is already in point_contacts and B is evaluated.
                        //If B is true, then CNT is not in point_contacts[Pa], thus the statements inside
                        //the if-statement are evaluated and Pa is added to point_contacts and CNTb is added
                        //to point_contacts_dist[Pa].
                        //
                        //If A and B are false, then Pa is already in point_contacts and CNTb is already
                        //in point_contacts[Pa] and C is evaluated. In this case, C is true only if the new
                        //junction distance is smaller than the one in point_contacts_dist[Pa][CNTb].
                        //Thus, it is only updated when a smaller junction is found
                        //
                        if (point_contacts.find(Pa) == point_contacts.end() || point_contacts[Pa].find(CNTb) == point_contacts[Pa].end() || point_contacts_dist[Pa][CNTb] > dist_junc) {
                            //hout<<"A="<<to_string(point_contacts.find(Pa) == point_contacts.end())<<" B="<<to_string(point_contacts[Pa].find(CNTb) == point_contacts[Pa].end())<<" C="<<to_string(point_contacts_dist[Pa][CNTb] > dist_junc)<<endl;
                            
                            //Update the junction distance
                            point_contacts_dist[Pa][CNTb] = dist_junc;
                            //hout<<"point_contacts_dist[Pa="<<Pa<<"][CNTb="<<CNTb<<"]="<< point_contacts_dist[Pa][CNTb]<<" Pb="<<Pb<<endl;
                            
                            //Add the correspoding point contacts
                            point_contacts[Pa][CNTb] = Pb;
                            //point_contacts[P2][CNT1] = P1;
                        }
                        //hout<<"point_contacts added"<<endl;
                    }
                }
            }
        }
    }
    
    return 1;
}
//This function merges labels in the labels vectors, so that labels_labels contains only proper labels
int Hoshen_Kopelman::Cleanup_labels(vector<int> &labels_labels, vector<int> &labels, int &n_labels)
{
    //Map of labels to know which ones are root labels and renumber them to a proper label
    vector<int> label_map(labels_labels.size(),-1);
    
    //This variable will be used to assign the proper label
    //It is initalized with 0, which is the first proper label
    int proper_label = 0;
    
    //First scan all the labels of labels to find the root labels
    //Create a map from root label to proper label
    //Proper labels are consecutive labels from 0 to N-1, where N is the number of clusters
    //i.e., these are merged and renumbered labels
    for (int i = 0; i < (int)labels_labels.size(); i++) {
        
        //if labels_labels_tmp[i] > 0, then i is a proper label
        if (labels_labels[i] > 0) {
            
            //Set label_label i with proper label counter
            label_map[i] = proper_label;
            
            //Increse the number of proper labels
            proper_label++;
        }
    }
    
    //Update the total number of proper labels
    n_labels = proper_label;
    
    //Scan all labels and change labels to keep only the proper labels
    for (int i = 0; i < (int)labels.size(); i++) {
        
        //Check that labels[i] is a valid label
        if (labels[i] != -1) {
            
            //Find the root label of labels[i]
            int root = Find_root(labels[i], labels_labels);
            
            //Find the proper label
            int proper = label_map[root];
            if (proper == -1) {
                hout << "Invalid proper label " << proper << ". Valid labels are positive integers or 0"<<endl;
                return 0;
            }
            
            //Change the current label by the proper label
            labels[i] = proper;
        }
    }
    
    return 1;
}
//This function renumbre GNP labels to have consecutive numbers after the last CNT label
int Hoshen_Kopelman::Renumber_gnp_labels(const int &cnt_labels, vector<int> &labels_gnp)
{
    //Scan all GNP labels
    for (size_t i = 0; i < labels_gnp.size(); i++) {
        
        //Check if label i is a proper label, i.e., it is positive or zero
        if (labels_gnp[i] >= 0) {
            
            //Add the number of CNT labels
            labels_gnp[i] = labels_gnp[i] + cnt_labels;
        }
    }
    
    return 1;
}
//This function compresses the "contact segments" for CNT-CNT contacts
int Hoshen_Kopelman::Compress_cnt_cnt_contact_segments(const Cutoff_dist &cutoffs, const map<long int, map<int, long int> > &point_contacts, const map<long int, map<int, double> > &point_contacts_dist, const vector<map<int, set<long int> > > &contact_elements)
{
    //Initialize the vector of elements
    elements_cnt.assign(contact_elements.size(), set<long int>());
    
    //Iterate over all CNTs
    for (int i = 0; i < (int)contact_elements.size(); i++) {
        
        //Check if CNT i has contacts with other CNTs
        if (!contact_elements[i].empty()) {
            
            //Iterate over the CNT in contact with CNT i
            for (map<int, set<long int> >::const_iterator i_map = contact_elements[i].begin(); i_map != contact_elements[i].end(); i_map++) {
                
                //Get the CNT number in contact with CNT i
                int CNTj = i_map->first;
                //hout<<"CNTi="<<i<<" CNTj="<<CNTj<<endl;
                
                //Start the iterator at the beginning of the set
                set<long int>::const_iterator i_set = i_map->second.begin();
                
                //Get the initial points in contact
                long int Pi0 = *i_set;
                //hout<<"Pi0="<<Pi0;
                long int Pj0 = point_contacts.at(Pi0).at(CNTj);
                //hout<<" Pj0="<<Pj0<<endl;
                
                //Get the initial junction distance
                double d_junction = point_contacts_dist.at(Pi0).at(CNTj);
                //hout<<"point_contacts_dist.at(Pi0="<<Pi0<<").at(CNTj="<<CNTj<<")="<<point_contacts_dist.at(Pi0).at(CNTj)<<" d_junction="<<d_junction<<endl;
                
                //Variables to store the "compress" contact segment
                long int Pi_junc = Pi0, Pj_junc = Pj0;
                double d_junc_min = d_junction;
                
                //Iterate over the points in CNT i in contact with CNTj and compress segements
                for (i_set++; i_set != i_map->second.end(); i_set++) {
                    
                    //Get the point in CNT i
                    long int Pi1 = *i_set;
                    
                    //Get the point in CNTj
                    //hout<<"Pi1="<<Pi1;
                    long int Pj1 = point_contacts.at(Pi1).at(CNTj);
                    //hout<<" Pj1="<<Pj1<<endl;
                    
                    //Get the junction distance
                    d_junction = point_contacts_dist.at(Pi1).at(CNTj);
                    //hout<<"d_junction="<<d_junction<<endl;
                    
                    //Check if the points are close enough to be in the same segment
                    if (Pi1 - Pi0 >= cutoffs.min_points || abs(Pj1 - Pj0) >= cutoffs.min_points) {
                        
                        //A new segment needs to be started
                        //But first, "compress" the previous segment
                        //That is, just add the shortest junction found to the
                        //vector of junctions
                        Junction j(Pi_junc, i, "CNT", Pj_junc, CNTj, "CNT", d_junc_min);
                        junctions_cnt.push_back(j);
                        //hout<<"compressed junction added in loop, Pi_junc="<<Pi_junc<<", CNTi="<<i<<", Pj_junc="<<Pj_junc<<", CNTj="<<CNTj<<" d_junc_min="<<d_junc_min<<endl;
                        
                        //Add the junction points to the vectors of elements
                        //hout<<"elements_cnt[i].size="<<elements_cnt[i].size()<<endl;
                        elements_cnt[i].insert(Pi_junc);
                        //hout<<"elements_cnt[CNTj].size="<<elements_cnt[CNTj].size()<<endl;
                        elements_cnt[CNTj].insert(Pj_junc);
                        
                        //Start the new segment by resetting the variables to the values
                        //of the junction that starts the new segment
                        d_junc_min = d_junction;
                        Pi_junc = Pi1;
                        Pj_junc = Pj1;
                    }
                    else {
                        
                        //Current contact is part of the same segment
                        //Check if junction distance and points in contact need updating
                        //Keep the points with the shortest junction distance since this
                        //is the path of lowest electrical resistance
                        if (d_junction < d_junc_min) {
                            
                            //A shorter distance was found, so update the variables for
                            //the junction for the compressed contact segment
                            d_junc_min = d_junction;
                            Pi_junc = Pi1;
                            Pj_junc = Pj1;
                        }
                    }
                    
                    //Update the value of Pi0 and Pj0
                    Pi0 = Pi1;
                    Pj0 = Pj1;
                }
                
                //"Compress" the last (or only) segment
                //That is, just add the shortest junction found in the last (or only) segment
                //to the vector of junctions
                Junction j(Pi_junc, i, "CNT", Pj_junc, CNTj, "CNT", d_junc_min);
                junctions_cnt.push_back(j);
                //hout<<"compressed junction added, Pi_junc="<<Pi_junc<<", CNTi="<<i<<", Pj_junc="<<Pj_junc<<", CNTj="<<CNTj<<" d_junc_min="<<d_junc_min<<endl;
                
                //Add the junction points to the vectors of elements
                //hout<<"elements_cnt[i].size="<<elements_cnt[i].size()<<endl;
                elements_cnt[i].insert(Pi_junc);
                //hout<<"elements_cnt[CNTj].size="<<elements_cnt[CNTj].size()<<endl;
                elements_cnt[CNTj].insert(Pj_junc);
            }
        }
    }
    
    return 1;
}
//This functions adds the first and last points of a CNT to the vector of elements
//These points were never added unless there is a junction on these points
int Hoshen_Kopelman::Complete_cnt_elements(const vector<vector<long int> > &structure_cnt)
{
    //Scan the vector of elements
    for (int i = 0; i < (int)elements_cnt.size(); i++) {
        
        //Check if element i is empty
        if (!elements_cnt[i].empty()) {
            
            //If the element is non-empty, then insert the first and last points
            elements_cnt[i].insert(elements_cnt[i].begin(), structure_cnt[i][0]);
            elements_cnt[i].insert(elements_cnt[i].end(), structure_cnt[i].back());
        }
    }
    
    return 1;
}
//This function makes the GNP clusters
int Hoshen_Kopelman::Make_gnp_clusters(const cuboid& sample, const vector<int> &gnps_inside, const vector<vector<int> > &sectioned_domain_gnp, const vector<GNP> &gnps, const Cutoff_dist& cutoffs, const int &n_labels_cnt, int &n_total_labels, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp, vector<int> &labels_gnp)
{
    //Vector of labels of labels
    vector<int> labels_labels_gnp;
    
    //Label the GNPs
    //hout<<"Label_gnps_in_window"<<endl;
    if (!Label_gnps_in_window(sample, gnps_inside, sectioned_domain_gnp, gnps, cutoffs, structure_gnp, points_gnp, labels_gnp, labels_labels_gnp)) {
        hout<<"Error in Make_gnp_clusters when calling Label_gnps_in_window"<<endl;
        return 0;
    }
    
    //Clean up the labels to find the proper labels, i.e. merged and consecutive labels starting at 0
    //hout<<"Cleanup_labels"<<endl;
    int n_gnp_labels = 0;
    if (!Cleanup_labels(labels_labels_gnp, labels_gnp, n_gnp_labels)) {
        hout << "Error in Make_gnp_clusters when calling Cleanup_labels" << endl;
        return 0;
    }
    
    //Check if there are CNT clusters
    if (n_labels_cnt) {
        
        //There are CNT clusters, thus rennumber the proper labels to have consecutive numbers
        //starting from n_total_labels
        if (!Renumber_gnp_labels(n_labels_cnt, labels_gnp)) {
            hout << "Error in Make_gnp_clusters when calling Renumber_gnp_labels" << endl;
            return 0;
        }
    }
    
    //Update the total number of labels
    n_total_labels = n_labels_cnt + n_gnp_labels;
    
    return 1;
}
//This function labels the GNPs so that they can be grouped into clusters
int Hoshen_Kopelman::Label_gnps_in_window(const cuboid& sample, const vector<int> &gnps_inside, const vector<vector<int> > &sectioned_domain_gnp, const vector<GNP> &gnps, const Cutoff_dist& cutoffs, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp, vector<int> &labels_gnp, vector<int> &labels_labels_gnp)
{
    //new_label will take the value of the newest cluster
    //Since GNP clusters are made after CNT clusters, the first label for GNPs is equal
    //to the number of CNT labels
    int new_label = 0;
    
    //Vector to check if a GNP-GNP contact has already been visited
    vector<set<int> > visited(gnps.size());
    
    //Scan every overlapping sub-region
    //hout<<"sectioned_domain_gnp.size()="<<sectioned_domain_gnp.size()<<endl;
    for (int i = 0; i < (int)sectioned_domain_gnp.size(); i++) 
    {
        //Size of subregion i
        int inner = (int)sectioned_domain_gnp[i].size();
        //hout<<"sectioned_domain_gnp[i="<<i<<"].size()="<<sectioned_domain_gnp[i].size()<<endl;
        
        //Scan all GNPs in subregion i
        for (int j = 0; j < inner - 1; j++) 
        {
            //Get GNP1
            int GNP1 = sectioned_domain_gnp[i][j];

            //Campare GNP1 with all other GNPs in sectioned domain i, if any
            for (int k = j + 1; k < inner; k++)
            {
                //Get GNP2
                int GNP2 = sectioned_domain_gnp[i][k];
                
                //Sort the GNP numbers
                int GNPa = min(GNP1, GNP2);
                int GNPb = max(GNP2, GNP2);
                //hout << "i=" << i << " j=" << j << " k=" << k << endl;
                //hout<<"GNPa="<<GNPa<<" GNPb="<<GNPb<<" GNP1="<<GNP1<<" GNP2="<<GNP2<<endl;
                
                //Check if the contact has already been visited
                if (visited[GNPa].find(GNPb) == visited[GNPa].end()) 
                {
                    //Contact GNPa-GNPb has not been visited

                    //Create a temporary GNP that will be moved in case there is interpenetration
                    GNP gnp_B = gnps[GNPb];

                    //Displacement variable in case there is interpenetration
                    //This is the displacement applied to GNPb in case there is interpenetration
                    Point_3D disp_tot;

                    //Penetration flag initialized with false
                    bool p_flag = false;

                    //Variables to store the distance between GNPa and GNPb and the direction 
                    //vect between them
                    double dist;
                    Point_3D N;

                    //Get the distance between GNPs and, in case there are interpenetrations,
                    //also get gnps[GNPb] in the new variable gnp_B
                    if (!Get_distance_between_GNPs(cutoffs, gnps, GNPa, GNPb, p_flag, disp_tot, N, dist, gnp_B))
                    {
                        hout << "Error in Label_gnps_in_window when calling Get_distance_between_GNPs" << endl;
                            return 0;
                    }

                    //Check if the separation between GNPs is below the cutoff for tunneling
                    //hout<<"GNPa="<<GNPa<<" GNPb="<<GNPb<<" dist="<<dist<<endl;
                    //hout<<"dist <= tunnel_cutoff = "<<(dist <= tunnel_cutoff)<<endl;
                    if (dist <= cutoffs.tunneling_dist) 
                    {
                        /* /Export GNPs that are found to be in contact for manually testing percolation
                        hout << "Contact: GNPa=" << GNPa << " GNPb=" << GNPb << " dist=" << dist << endl;
                        VTK_Export VTK;
                        VTK.Export_single_gnp(gnps[GNPa], "gnp_" + to_string(GNPa) + ".vtk");
                        VTK.Export_single_gnp(gnps[GNPb], "gnp_" + to_string(GNPb) + ".vtk");// */

                        //Junction points on GNPs
                        Point_3D PointA, PointB;
                        
                        //Find the points of contact on each GNP, and add them if both
                        //are inside the sample
                        //hout<<"Add_junction_points_for_gnps"<<endl;
                        if (!Find_junction_points_for_gnps(sample, gnps[GNPa], gnp_B, N, dist, points_gnp, PointA, PointB)) 
                        {
                            hout << "Error in Label_gnps_in_window when calling Add_junction_points_for_gnps" << endl;
                            return 0;
                        }

                        //Check if there is interpenetration of GNPs
                        if (p_flag)
                        {
                            //There is interpenetration of GNPs, so GNPb was moved
                            //Thus, move the point in GNPb towards GNPb's original position
                            PointB = PointB - disp_tot;
                            //Flag might be lost in the operation above
                            PointB.flag = GNPb;
                        }

                        //Check if both junction points are inside the sample
                        if (!PointA.is_outside_cuboid(sample) && !PointB.is_outside_cuboid(sample))
                        {
                            //Junction is inside the sample so assign the GNPs into the same cluster
                            //Here is where the actual HK76 algorithm takes place
                            if (!HK76(GNPa, GNPb, new_label, labels_gnp, labels_labels_gnp)) 
                            {
                                hout << "Error in Label_gnps_in_window when calling HK76" << endl;
                                return 0;
                            }

                            //Add junction to vector of junctions
                            if (!Add_gnp_junction(gnps, GNPa, GNPb, p_flag, PointA, PointB, dist, cutoffs.tol_gnp2, points_gnp, structure_gnp))
                            {
                                hout << "Error in Label_gnps_in_window when calling Add_gnp_junction" << endl;
                                return 0;
                            }
                        }
                        //else {
                            //The junction is actually ouside the sample, so it should be ignored
                            //hout << "Ignore junction between GNP " << GNP_A.flag << " and GNP " << GNP_B.flag << endl;
                        //}
                    }
                    
                    //Mark the contact as visited
                    visited[GNPa].insert(GNPb);
                }
            }
        }
    }
    
    return 1;
}
//This function moves GNPb to a position where it does not penetrate GNPa, however the GNP in its
//new location is saved into a new variable
int Hoshen_Kopelman::Get_distance_between_GNPs(const Cutoff_dist& cutoffs, const vector<GNP>& gnps, const int& GNPa, const int& GNPb, bool& p_flag, Point_3D& disp_tot, Point_3D& N, double& dist, GNP& gnp_B)
{
    //Collision detection object
    Collision_detection CL;

    //Variables for the GJK
    vector<Point_3D> simplex;

    //Use the GJK for distance to calculate the separation between GNPs
    //and the direction vector
    //The distance between GNPs is saved in the output variable dist
    //hout<<"CL.GJK_distance"<<endl;
    if (!CL.GJK_distance(gnps[GNPa], gnps[GNPb], simplex, dist, N, p_flag)) {
        hout << "Error in Get_distance_between_GNPs when calling CL.GJK_distance" << endl;
        return 0;
    }

    //Check if the GNPs are interpenetrating each other
    if (p_flag)
    {
        //GNPs are interpenetrating each other
        //Use EPA to find the penetration depth PD and direction vector N
        if (!CL.EPA(gnps[GNPa].vertices, gnp_B.vertices, simplex, N, dist)) {
            hout << "Error in Get_distance_between_GNPs when calling EPA" << endl;
            string gnp1_str = "gnp_" + to_string(gnps[GNPa].flag) + ".vtk";
            string gnp2_str = "gnp_" + to_string(gnps[GNPb].flag) + ".vtk";
            hout << "GNPs are exported into files " << gnp1_str << " and " << gnp2_str << endl;
            //Export GNPs thta caused the error
            VTK_Export VTK;
            VTK.Export_single_gnp(gnps[GNPa], gnp1_str);
            VTK.Export_single_gnp(gnps[GNPb], gnp2_str);
            return 0;
        }

        //Calculate total displacement for GNPb
        disp_tot = N * (dist + cutoffs.van_der_Waals_dist);

        //Move gnp_B so that there is no interpenetration and so that the 
        //GNPs are separated the van der Waals distance
        //Use a Generate_Network object
        Generate_Network GN;
        if (!GN.Move_gnp(disp_tot, gnp_B))
        {
            hout << "Error in Get_distance_between_GNPs when calling GN.Move_gnp" << endl;
            return 0;
        }

        //Now GNPs are separated a distance equal to the van der Waals distance
        //So set the distance to that value
        dist = cutoffs.van_der_Waals_dist;
    }

    return 1;
}
//This function finds the points in two GNPs that will serve as junction points
//These points are also used to define the GNP resistors
int Hoshen_Kopelman::Find_junction_points_for_gnps(const cuboid &sample, const GNP &GNP_A, const GNP &GNP_B, const Point_3D &N, const double &distance, vector<Point_3D> &points_gnp, Point_3D& PointA, Point_3D& PointB)
{
    //Vectors to save the simplices in A and B that are closest to each other
    vector<int> simplexA, simplexB;
    int v_sumA = 0, v_sumB = 0;
    if (!Find_closest_simplices_of_gnps_in_contact(GNP_A, GNP_B, N, distance, simplexA, simplexB, v_sumA, v_sumB)) {
        hout<<"Error in Add_junction_points_for_gnps when calling Find_closest_simplices_of_gnps_in_contact"<<endl;
        return 0;
    }
    //hout<<"simplexA.size="<<simplexA.size()<<" simplexB.size="<<simplexB.size()<<" v_sumA="<<v_sumA<<" v_sumB="<<v_sumB<<endl;
    
    //Using the simplices in GNP_A and GNP_B find the two points that define a junction between
    //the two GNPs
    if (simplexA.size() <= simplexB.size()) {
        
        //Simplex A is the one with less vertices or it has the same number of vertices as B
        if (!Find_junction_points_in_gnps(simplexA, simplexB, GNP_A, GNP_B, N, distance, v_sumA, v_sumB, PointA, PointB)) {
            hout<<"Error in Add_junction_points_for_gnps when calling Find_junction_points_in_gnps (1)"<<endl;
            return 0;
        }
    }
    else {
        
        //Simplex B is the one with less vertices
        //Invert the direction of N
        //Interchange GNP_A and GNP_B
        if (!Find_junction_points_in_gnps(simplexB, simplexA, GNP_B, GNP_A, N*(-1.0), distance, v_sumB, v_sumA, PointB, PointA)) {
            hout<<"Error in Add_junction_points_for_gnps when calling Find_junction_points_in_gnps (2)"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function finds the two simplices that are closest to each other when two GNPs are in
//electrical contact
int Hoshen_Kopelman::Find_closest_simplices_of_gnps_in_contact(const GNP &GNP_A, const GNP &GNP_B, const Point_3D &N, const double &distance, vector<int> &simplexA, vector<int> &simplexB, int &v_sumA, int &v_sumB)
{
    //Sets to store the vertices of the simplex in the Minkowski sum
    //Sets are used to avoid repetition of vertices
    set<int> setA, setB;
    
    //Vectors to store points that help me debug the code
    //vector<Point_3D> cnt_tmp, cnt_simplex;
    
    //Iterate over all vertices of the Minkowski sum
    //i iterates over the vertices of A
    for (int i = 0; i < 8; i++) {
        //j iterates over the vertices of B
        for (int j = 0; j < 8; j++) {
            
            //Calculate the point in the Minkowski sum
            Point_3D Q = GNP_A.vertices[i] - GNP_B.vertices[j];
            //cnt_tmp.push_back(Q);
            //hout<<"Q_"<<i<<j<<"="<<Q.str(10)<<" A="<<GNP_A.vertices[i].str(10)<<" B="<<GNP_B.vertices[j].str(10)<<endl;
            
            //Check if Q is part of the simplex closest to the origin by checking if the
            //cutoff is computationally zero, i.e., its absolute value is less than Zero
            //double tmp = abs(Q.dot(N));
            //hout<<"abs(Q.dot(N))["<<tmp<<"] - distance["<<distance<<"]="<<abs(tmp - distance)<<endl;
            if (abs(abs(Q.dot(N)) - distance) <= Zero) {
                
                //Q is at the same distance to the origin as the simplex in the Minkowski sum
                //that is closest to the origin
                //Thus, vertex i in GNP_A and vertex j in GNP_B are part of the simplices
                //In those GNPs that are closest to each other
                setA.insert(i);
                setB.insert(j);
                //double tmp = abs(Q.dot(N));
                //hout<<"abs(Q.dot(N))["<<tmp<<"] - distance["<<distance<<"]="<<abs(tmp - distance)<<endl;
                //hout<<" i="<<i<<"->simplexA j="<<j<<"->simplexB"<<endl<<endl;
                //cnt_simplex.push_back(Q);
            }
        }
    }
    
    /* /hout<<"=========================="<<endl<<endl;
    hout<<endl<<"GNPs:"<<endl;
    //Iterates over the vertices of A
    for (int i = 0; i < 8; i++) {
        hout<<GNP_A.vertices[i].str(10)<<endl;
    }
    //Iterates over the vertices of B
    for (int i = 0; i < 8; i++) {
        hout<<GNP_B.vertices[i].str(10)<<endl;
    }
    hout<<endl<<endl;*/
    
    /*VTK_Export VTK;
    vector<GNP> tmp_vec;
    GNP tmp_gnp;
    //i iterates over the vertices of B
    for (int i = 0; i < 8; i++) {
        //j iterates over the vertices of A
        for (int j = 0; j < 8; j++) {
            //Create a new GNP
            tmp_gnp.vertices[j] = GNP_A.vertices[j] - GNP_B.vertices[i];
        }
        tmp_vec.push_back(tmp_gnp);
    }*/
    /*VTK_Export VTK;
    vector<GNP> tmp_vec(2);
    tmp_vec[0] = GNP_A; tmp_vec[1] = GNP_B;
    VTK.Export_gnps(tmp_vec, "gnps_debug.vtk");*/
    //
    //
    //VTK.Export_single_cnt(cnt_tmp, "cnt_tmp.vtk");
    //VTK.Export_single_cnt(cnt_simplex, "cnt_simplex.vtk");
    
    //Fill the vectors of the simplices using the sets
    simplexA.assign(setA.begin(), setA.end());
    simplexB.assign(setB.begin(), setB.end());
    
    //Calculate the vertices sums
    for (size_t i = 0; i < simplexA.size(); i++) {
        v_sumA = v_sumA + max(1,simplexA[i]);
    }
    for (size_t i = 0; i < simplexB.size(); i++) {
        v_sumB = v_sumB + max(1,simplexB[i]);
    }
    
    return 1;
}
//This function finds the two points that define the junction between two GNPs
//It is assumed that simplexA has less or the same number of vertices as simplexB
//It is also assumed that the direction vector goes from GNP_A to GNP_B
int Hoshen_Kopelman::Find_junction_points_in_gnps(const vector<int> &simplexA, const vector<int> &simplexB, const GNP &GNP_A, const GNP &GNP_B, const Point_3D &N, const double &distance, const int &face_sumA, const int &face_sumB, Point_3D &PointA, Point_3D &PointB)
{
    //Find the case of determining the junction points
    if (simplexA.size() == 1) {
        
        //Vertex in simplexA
        int vA = simplexA[0];
        
        //PointA is the vertex in simplexA
        PointA = GNP_A.vertices[vA];
        
        //Find PointB, depending on the size of simplexB
        if (!Find_point_b_for_vertex_in_simplex_a(simplexB, GNP_B, N, distance, PointA, PointB)) {
            hout << "GNP_A=" << GNP_A.flag << " GNP_B=" << GNP_B.flag << endl;
            hout<<"Error in Find_junction_points_in_gnps when calling Find_point_b_for_vertex_in_simplex_a"<<endl;

            //GNP_A could not be saved within fucntion Find_point_b_for_vertex_in_simplex_a since
            //it is not one of it arguments
            //Thus it is saved here
            VTK_Export VTK;
            VTK.Export_single_gnp(GNP_A, "gnp_A.vtk");

            return 0;
        }
    }
    else if (simplexA.size() == 2) {
        
        //Find PointB, depending on the size of simplexB
        if (!Find_point_b_for_edge_in_simplex_a(simplexA, simplexB, face_sumB, GNP_A, GNP_B, N, distance, PointA, PointB)) {
            hout<<"Error in Find_junction_points_in_gnps when calling Find_point_b_for_edge_in_simplex_a"<<endl;
            return 0;
        }
    }
    else if (simplexA.size() == 4) {
        
        //Find PointB, depending on the size of simplexB
        if (!Find_point_b_for_face_in_simplex_a(simplexA, simplexB, face_sumA, face_sumB, GNP_A, GNP_B, N, distance, PointA, PointB)) {
            hout<<"Error in Find_junction_points_in_gnps when calling Find_point_b_for_face_in_simplex_a"<<endl;
            return 0;
        }
    }
    else {
        hout<<"Error in Find_junction_points_in_gnps: simplexA has an invalid size="<<simplexA.size()<<", it size must be 1, 2, or 4"<<endl;
        hout<<"GNP for simplex of size " << simplexA.size() << " is saved in file gnp_A.vtk." << endl;
        hout<<"Second GNP is saved in file gnp_B.vtk." << endl;
        hout << "GNP_A=" << GNP_A.flag << " GNP_B=" << GNP_B.flag << endl;

        VTK_Export VTK;
        VTK.Export_single_gnp(GNP_A, "gnp_A.vtk");
        VTK.Export_single_gnp(GNP_B, "gnp_B.vtk");

        return 0;
    }
    
    //Set the flag of PointA
    PointA.flag = GNP_A.flag;
    
    //Set the flag of PointB
    PointB.flag = GNP_B.flag;
    
    return 1;
}
//This function finds the contact points when there is a vertex-vertex, vertex-edge or vertex-face contact
int Hoshen_Kopelman::Find_point_b_for_vertex_in_simplex_a(const vector<int> &simplexB, const GNP &GNP_B, const Point_3D &N, const double &distance, const Point_3D &PointA, Point_3D &PointB)
{
    
    //Find point B depending on the number of vertices in simplexB
    if (simplexB.size() == 1) {
        
        //Vertex in simplexB
        int vB = simplexB[0];
        
        //PointB is the vertex in simplexB
        PointB = GNP_B.vertices[vB];
    }
    else if (simplexB.size() == 2 || simplexB.size() == 4) {
        
        //To find the projection on the edge or face in simplexB, I just need to move PointA
        //towards the edge or face by the given distance and direction
        PointB = PointA + N*distance;
        
    }
    else {
        hout<<"Error in Find_point_b_for_vertex_in_simplex_a: simplexB has an invalid size="<<simplexB.size()<<", it size must be 1, 2, or 4"<<endl;
        hout << "GNP for simplex of size 1 is saved in file gnp_A.vtk." << endl;
        hout << "GNP for simplex of size " << simplexB.size() << " is saved in file gnp_B.vtk." << endl;

        //Export the GNPs for debugging
        VTK_Export VTK;
        VTK.Export_single_gnp(GNP_B, "gnp_B.vtk");
        return 0;
    }
    
    return 1;
}
//This function finds the contact points when there is an edge-edge or edge-face contact
int Hoshen_Kopelman::Find_point_b_for_edge_in_simplex_a(const vector<int> &simplexA, const vector<int> &simplexB, const int &face_sum, const GNP &GNP_A, const GNP &GNP_B, const Point_3D &N, const double &distance, Point_3D &PointA, Point_3D &PointB)
{
    
    //Calculate the displacement from A to B
    Point_3D disp = N*distance;
    
    //Move simplexA towards simplexB
    Point_3D edgeA[] = {
        GNP_A.vertices[simplexA[0]] + disp,
        GNP_A.vertices[simplexA[1]] + disp
    };
    
    //Find points A and B depending on the number of vertices in simplexB
    if (simplexB.size() == 2) {
        
        //Go to the case where both simplexA and simplexB are an edge
        if (!Find_point_for_edge_edge_contact(GNP_A, GNP_B, edgeA, simplexB, disp, PointA, PointB)) {
            hout<<"Error in Find_point_b_for_edge_in_simplex_a when calling Find_point_for_edge_edge_contact"<<endl;
            return 0;
        }
    }
    else if (simplexB.size() == 4) {
        
        //Go to the case where simplexB is a face
        if (!Find_point_for_edge_face_contact(face_sum, GNP_A, GNP_B, edgeA, disp, PointA, PointB)) {
            hout<<"Error in Find_point_b_for_edge_in_simplex_a when calling Find_point_for_edge_face_contact"<<endl;
            return 0;
        }
    }
    else {
        hout<<"Error in Find_point_b_for_edge_in_simplex_a: simplexB has an invalid size="<<simplexB.size()<<", it size must be 2 or 4."<<endl;
        hout << "GNP for simplex of size 2 is saved in file gnp_A.vtk." << endl;
        hout << "GNP for simplex of size "<< simplexB.size() <<" is saved in file gnp_B.vtk." << endl;
        hout << "GNP_A=" << GNP_A.flag << " GNP_B=" << GNP_B.flag << endl;

        //Export the GNPs for debugging
        VTK_Export VTK;
        VTK.Export_single_gnp(GNP_A, "gnp_A.vtk");
        VTK.Export_single_gnp(GNP_B, "gnp_B.vtk");

        return 0;
    }
    return 1;
}
//This function finds the contact point for an edge-edge contact, which can be:
//a) The intersection of the two edges
//b) The midpoint of the overlapping segement of the two edges
int Hoshen_Kopelman::Find_point_for_edge_edge_contact(const GNP &GNP_A, const GNP &GNP_B, const Point_3D edgeA[], const vector<int> &simplexB, const Point_3D &disp, Point_3D &PointA, Point_3D &PointB)
{
    //Get vectors in the direction of each edge
    Point_3D eA = edgeA[1] - edgeA[0];
    Point_3D eB = GNP_B.vertices[simplexB[1]] - GNP_B.vertices[simplexB[0]];
    
    //Check if edges are colinear
    Point_3D cr = eA.cross(eB);
    if (abs(cr.x) < Zero && abs(cr.y) < Zero && abs(cr.z) < Zero) {
        
        //Cross product of vectors along the edges is zero, thus edges are colinear
        
        //Find the midpoint of the intersection of the edges defined by the vertices in
        //simplices A and B
        
        //Lambda function to calculate the lambda value of a point
        auto calc_lambda = [](const double &x1, const double &x2, const double &x_new){
            return ((x_new - x1)/(x2 - x1));
        };
        
        //Calculate the lambdas of the vertices of the other edge
        double lambda0 = calc_lambda(edgeA[0].x, edgeA[1].x, GNP_B.vertices[simplexB[0]].x);
        double lambda1 = calc_lambda(edgeA[0].x, edgeA[1].x, GNP_B.vertices[simplexB[1]].x);
        
        //One of the lambdas will be in the range [0,1]
        //The other will be either greater than 1 or less than 0
        //Only if the two edges are perfectly overlappin, one lambda will be
        //0 and the other 1, so also check for that case
        if ( (abs(lambda0) < Zero && abs(lambda1 - 1) < Zero) ||
             (abs(lambda1) < Zero && abs(lambda0 - 1) < Zero) ) {
            
            //If the edges are perfectly overlappig, just take the average of the vertices
            PointB = (GNP_B.vertices[simplexB[0]] + GNP_B.vertices[simplexB[1]])/2;
        }
        if (lambda0 >= Zero && lambda0 - 1 <= Zero) {
            
            //Vertex 0 in simplexB is on the translated edgeA
            
            //Calculate point on translated edgeA
            Point_3D S0 = edgeA[0] + eA*lambda0;
            
            //Define PointB based on the sign of lambda 1
            PointB = (lambda1 < Zero)? (edgeA[0] + S0)/2 : (edgeA[1] + S0)/2;
        }
        else {
            
            //Vertex 1 in simplexB is on the translated edgeA
            
            //Calculate point on translated edgeA
            Point_3D S1 = edgeA[0] + eA*lambda1;
            
            //Define PointB based on the sign of lambda 1
            PointB = (lambda0 < Zero)? (edgeA[0] + S1)/2 : (edgeA[1] + S1)/2;
        }
    }
    else {
        
        //Edges are not colinear, so find the intersection of the edges
        
        //Find the lambda value of the intersection using edgeB (simplexB)
        //Lambda = 0 at GNP_B.vertices[simplexB[0]]
        //P = GNP_B.vertices[simplexB[0]] +
        //      (GNP_B.vertices[simplexB[1]] - GNP_B.vertices[simplexB[0]])*lambda
        double lambda = Lambda_of_two_lines(GNP_B.vertices[simplexB[1]], GNP_B.vertices[simplexB[0]], edgeA[0], edgeA[1]);
        
        //Calculate the intersection on edgeB (simplexB)
        PointB = GNP_B.vertices[simplexB[0]] + eB*lambda;
    }
    
    //Calculate PointA by using the calculated displacement
    PointA = PointB - disp;
    
    return 1;
}
//This function finds the contact point for an edge-face contact
int Hoshen_Kopelman::Find_point_for_edge_face_contact(const int &face_sum, const GNP &GNP_A, const GNP &GNP_B, const Point_3D edgeA[], const Point_3D &disp, Point_3D &PointA, Point_3D &PointB)
{
    //Find centroid of the intersection of the faces defined by the vertices in
    //simplices A and B
    
    //Get the edges and normals of simplexB
    vector<Edge> edges(4);
    vector<int> normals(4);
    if (!Get_edges_of_face(face_sum, edges, normals)) {
        hout<<"Error in Find_point_for_edge_face_contact when calling Get_edges_of_face"<<endl;
        return 0;
    }
    
    //Variables to count the negative dot products
    int countA0 = 0, countA1 = 0;
    
    //Arrays to store positive dot products
    int arrA0[] = {-1,-1};
    int arrA1[] = {-1,-1};
    
    //Go though all edges of faceB and find the one that intersects the edgeA
    for (int i = 0; i < 4; i++) {
        
        //Vertex 1 of current edge
        int v1i = edges[i].v1;
        
        //Current normal
        int Ni = normals[i];
        
        //Check in which side of the edge the two points are
        int dotA0 = (edgeA[0] - GNP_B.vertices[v1i]).dot(GNP_B.faces[Ni].N) < Zero;
        int dotA1 = (edgeA[1] - GNP_B.vertices[v1i]).dot(GNP_B.faces[Ni].N) < Zero;
        
        //Check if the count of negative dot products for vertex A0 needs to be updated
        if (dotA0) {
            //edgeA[0] is likely inside the face
            countA0++;
        }
        else {
            //Save the edge in B that has results in a positive dot product
            //There can be at most two positive dot prodects, so this if-statement is enough
            if (arrA0[0] == -1) {
                arrA0[0] = i;
            }
            else {
                arrA0[1] = i;
            }
        }
        //Check if the count of negative dot products for vertex A0 needs to be updated
        if (dotA1) {
            //edgeA[1] is likely inside the face
            countA1++;
        }
        else {
            //Save the edge in B that has results in a positive dot product
            //There can be at most two positive dot prodects, so this if-statement is enough
            if (arrA1[0] == -1) {
                arrA1[0] = i;
            }
            else {
                arrA1[1] = i;
            }
        }
    }
    
    //If edgeA is not entirely inside faceB (simplexB), then we need
    //to get the intersection points. If a vertex is inside faceB then
    //the intersection point is the same as the vertex
    //So, according to the count of negative dot products, get the intersection points
    
    //Get intersection point from vertex A0
    Point_3D A0 = Get_intersection_with_face_edge(countA0, arrA0, edgeA[0], edgeA[1], edges, GNP_B);
    
    //Check for an error
    if (A0.flag == -1) {
        hout<<"Error in Find_point_for_edge_face_contact when calculating an intersecting point (1)"<<endl;
        return 0;
    }
    
    //Get intersection point from vertex A1
    Point_3D A1 = Get_intersection_with_face_edge(countA1, arrA1, edgeA[1], edgeA[0], edges, GNP_B);
    
    //Check for an error
    if (A1.flag == -1) {
        hout<<"Error in Find_point_for_edge_face_contact when calculating an intersecting point (2)"<<endl;
        return 0;
    }
    
    //Get the average of the intersection points
    PointB = (A0 + A1)*0.5;
    
    //To calculate PointA, move PointB towards simplexA
    PointA = PointB - disp;
    return 1;
}
//This function generates the vector of edges and the vector of normals of a given GNP face
int Hoshen_Kopelman::Get_edges_of_face(const int &face_sum, vector<Edge> &edges, vector<int> &normals)
{
    //Make the edges vector depending on the face
    switch (face_sum) {
        case 7:
            edges[0].v1 = 0;
            edges[0].v2 = 3;
            edges[1].v1 = 0;
            edges[1].v2 = 1;
            edges[2].v1 = 1;
            edges[2].v2 = 2;
            edges[3].v1 = 2;
            edges[3].v2 = 3;
            normals[0] = 2;
            normals[1] = 3;
            normals[2] = 4;
            normals[3] = 5;
            break;
        case 22:
            edges[0].v1 = 4;
            edges[0].v2 = 7;
            edges[1].v1 = 4;
            edges[1].v2 = 5;
            edges[2].v1 = 5;
            edges[2].v2 = 6;
            edges[3].v1 = 6;
            edges[3].v2 = 7;
            normals[0] = 2;
            normals[1] = 3;
            normals[2] = 4;
            normals[3] = 5;
            break;
        case 15:
            edges[0].v1 = 0;
            edges[0].v2 = 3;
            edges[1].v1 = 3;
            edges[1].v2 = 7;
            edges[2].v1 = 7;
            edges[2].v2 = 4;
            edges[3].v1 = 0;
            edges[3].v2 = 4;
            normals[0] = 0;
            normals[1] = 5;
            normals[2] = 1;
            normals[3] = 3;
            break;
        case 11:
            edges[0].v1 = 1;
            edges[0].v2 = 0;
            edges[1].v1 = 0;
            edges[1].v2 = 4;
            edges[2].v1 = 4;
            edges[2].v2 = 5;
            edges[3].v1 = 5;
            edges[3].v2 = 1;
            normals[0] = 0;
            normals[1] = 2;
            normals[2] = 1;
            normals[3] = 4;
            break;
        case 14:
            edges[0].v1 = 1;
            edges[0].v2 = 2;
            edges[1].v1 = 2;
            edges[1].v2 = 6;
            edges[2].v1 = 6;
            edges[2].v2 = 5;
            edges[3].v1 = 5;
            edges[3].v2 = 1;
            normals[0] = 0;
            normals[1] = 5;
            normals[2] = 1;
            normals[3] = 3;
            break;
        case 18:
            edges[0].v1 = 2;
            edges[0].v2 = 3;
            edges[1].v1 = 3;
            edges[1].v2 = 7;
            edges[2].v1 = 7;
            edges[2].v2 = 6;
            edges[3].v1 = 6;
            edges[3].v2 = 2;
            normals[0] = 0;
            normals[1] = 2;
            normals[2] = 1;
            normals[3] = 4;
            break;
        
        default:
            hout<<"Error in Get_edges_of_face. Invalid face number. Face number goes from 0 to 5. Input was "<<face_sum<<endl;
            return 0;
    }
    
    return 1;
}
//This function determines the intersecting point of an edge from the side of a given vertex
//It is assumed A0 is this vertex in question
Point_3D Hoshen_Kopelman::Get_intersection_with_face_edge(const int &countA, const int arrA[], const Point_3D &A0, const Point_3D &A1, const vector<Edge> &edges, const GNP GNP_B)
{
    //Point to store the solution
    Point_3D P;
    //Set the flag to zero as it will be used to check for errors
    P.flag = 0;
    
    //Check the count of negative dot products
    if (countA == 4) {
        
        //Vertex of edgeA is inside faceB, thus return this vertex (i.e., vertex A0)
        return A0;
    }
    else if (countA  == 3) {
        
        //Vertex A0 is closest to an edge of faceB
        //This edge is stored in arrA[0], so get the edge number
        int eB = arrA[0];
        
        //Get the vertices of the edges
        int v1 = edges[eB].v1;
        int v2 = edges[eB].v2;
        
        //Calculate the lambda of the equation using the two vertices in edgeA
        double lambda = Lambda_of_two_lines(A0, A1, GNP_B.vertices[v1], GNP_B.vertices[v2]);
        
        //Check that lambda is in the range [0,1]
        if ( (lambda < Zero && abs(lambda) > Zero) || lambda > 1.0) {
            hout<<"Error in Get_intersection_with_face_edge (1): lambda is outside the range [0,1], lambda="<<lambda<<endl;
            P.flag = -1;
            return P;
        }
        
        //Use lambda to calculate the intersecting point
        P = A1 + (A0 - A1)*lambda;
    }
    else if (countA == 2) {
        
        //Vertex A0 is closest to a vertex of faceB
        //This vertex is the intersection of two edges and are stored in arrA[0] and arrA[1],
        //so get the edge numbers
        int eB = arrA[0];
        
        //Get the vertices of the edges
        int v1 = edges[eB].v1;
        int v2 = edges[eB].v2;
        
        //Calculate the lambda of the equation using the two vertices in edgeA
        double lambda = Lambda_of_two_lines(A0, A1, GNP_B.vertices[v1], GNP_B.vertices[v2]);
        
        //Check the sign of lambda
        if (lambda < Zero && abs(lambda) > Zero) {
            
            //lambda is negative, so edgeA intersects the straight line that contains
            //edge eB1, but this intersecting point is not contained in edge eB1
            //Thus, get the second edge shared by the vertex of faceB closest to A0
            eB = arrA[1];
            
            //Get the vertices of the edges
            v1 = edges[eB].v1;
            v2 = edges[eB].v2;
            
            //Calculate the lambda of the equation using the two vertices in edgeA
            lambda = Lambda_of_two_lines(A0, A1, GNP_B.vertices[v1], GNP_B.vertices[v2]);
            
            //Check that lambda is in the range [0,1]
            //This time it cannot be outside this range, otherwise it means edgeA
            //actually does not intersect any of these edges in faceB
            if ( (lambda < Zero && abs(lambda) > Zero) || lambda > 1.0) {
                hout<<"Error in Get_intersection_with_face_edge (2): lambda is outside the range [0,1], lambda="<<lambda<<endl;
                P.flag = -1;
                return P;
            }
        }
        
        //Use the positive lambda to calculate the intersecting point
        P = A1 + (A0 - A1)*lambda;
        
    }
    else {
        hout<<"Error in Get_intersection_with_face_edge: invalid value of countA. It can only be 2, 3, or 4. countA="<<countA<<endl;
        P.flag = -1;
        return P;
    }
    
    return P;
}
//This function calcualtes the lambda value of the intersection of two lines defined by two points
//The lambda value corresponds to the line equation defined by the two first points in
//arguments of the function
//Lambda = 0 at P_out
double Hoshen_Kopelman::Lambda_of_two_lines(const Point_3D &P_in, const Point_3D &P_out, const Point_3D &edge1, const Point_3D &edge2)
{
    //Calculate det(A)
    double detA = (edge2.x - edge1.x)*(P_in.y - P_out.y)
    - (P_in.x - P_out.x)*(edge2.y - edge1.y);
    
    //Calcualte det(A_lambda_P)
    double detAP = (edge2.x - edge1.x)*(edge1.y - P_out.y)
    - (edge1.x - P_out.x)*(edge2.y - edge1.y);
    
    //Calculate the lambda value
    return detAP/detA;
}
//This function finds the contact points when there is a face-face contact
int Hoshen_Kopelman::Find_point_b_for_face_in_simplex_a(const vector<int> &simplexA, const vector<int> &simplexB, const int &face_sumA, const int &face_sumB, const GNP &GNP_A, const GNP &GNP_B, const Point_3D &N, const double &distance, Point_3D &PointA, Point_3D &PointB)
{
    //Check the size of simplexB
    if (simplexB.size() == 4) {
        
        //Find the projection of PointA in the surface defined by the four vertices in simplexB
        
        //Get the vertices and normals from simplexB
        vector<Edge> edgesB(4);
        vector<int> normalsB(4);
        if (!Get_edges_of_face(face_sumB, edgesB, normalsB)) {
            hout<<"Error in Find_point_b_for_edge_in_simplex_a when calling Get_edges_of_face B"<<endl;
            return 0;
        }
        
        //Get the vertices in A that are inside or in the boundary of B
        vector<int> verticesA_inside;
        if (!Get_vertices_inside_face(simplexA, edgesB, normalsB, GNP_A, GNP_B, verticesA_inside)) {
            hout<<"Error in Find_point_b_for_edge_in_simplex_a when calling Get_vertices_inside_face B"<<endl;
            return 0;
        }
        
        //Check if the four vertices of A were added to the vector of vertices
        if (verticesA_inside.size() == 4) {
            
            //Faces of GNP_A and GNP_B are perfectly aligned and oriented
            
            //Set point A and B to zero
            PointA.set(0, 0, 0);
            PointB.set(0, 0, 0);
            
            //Get the average point of each GNP face in contact
            for (int k = 0; k < (int)simplexB.size(); k++) {
                
                //Get vertex in A and add it to PointA
                int vA = simplexA[k];
                PointA = PointA + GNP_A.vertices[vA];
                
                //Get vertex in B and add it to PointB
                int vB = simplexB[k];
                PointB = PointB + GNP_B.vertices[vB];
            }
            
            //Calcualte the averages
            PointA = PointA/4;
            PointB = PointB/4;
        }
        //If not all vertices of A are on face B, then get the vertices of B on face A
        else {
            
            //Get the vertices and normals from simplexA
            vector<Edge> edgesA(4);
            vector<int> normalsA(4);
            if (!Get_edges_of_face(face_sumA, edgesA, normalsA)) {
                hout<<"Error in Find_point_b_for_edge_in_simplex_a when calling Get_edges_of_face A"<<endl;
                return 0;
            }
            
            //Get the vertices in B that are inside or in the boundary of A
            vector<int> verticesB_inside;
            if (!Get_vertices_inside_face(simplexB, edgesA, normalsA, GNP_B, GNP_A, verticesB_inside)) {
                hout<<"Error in Find_point_b_for_edge_in_simplex_a when calling Get_vertices_inside_face A"<<endl;
                return 0;
            }
            
            //Calcualte the displacement
            Point_3D disp = N*distance;
            
            //Check the case when only one vertex is inside the other face
            if (verticesA_inside.size() == 0 && verticesB_inside.size() == 1) {
                
                //Initialize PointB with the vertes inside face A
                PointB = GNP_B.vertices[verticesB_inside[0]];
                
                //Find intersections with edge in A and get average PointA
                if (!Average_point_with_two_interserctions(edgesA, normalsA, GNP_A, edgesB, normalsB, GNP_B, disp, verticesB_inside[0], PointB)) {
                    hout<<"Error in Find_point_b_for_edge_in_simplex_a when calling Average_point_with_two_interserctions (1)"<<endl;
                    return 0;
                }
                
                //Move B towards A to obtain the point in A
                PointA = PointB - disp;
            }
            else if (verticesA_inside.size() == 1 && verticesB_inside.size() == 0) {
                
                //Initialize PointA with the vertes inside face B
                PointA = GNP_A.vertices[verticesA_inside[0]];
                
                //Find intersections with edge in B and get average PointB
                if (!Average_point_with_two_interserctions(edgesB, normalsB, GNP_B, edgesB, normalsB, GNP_B, disp*(-1.0), verticesA_inside[0], PointA)) {
                    hout<<"Error in Find_point_b_for_edge_in_simplex_a when calling Average_point_with_two_interserctions (2)"<<endl;
                    return 0;
                }
                
                //Move A towards B to obtain the point in B
                PointB = PointA + disp;
            }
            else {
                
                //Set B to zero
                PointB.set(0, 0, 0);
                
                //Add displaced vertices of A inside B
                for (int k = 0; k < (int)verticesA_inside.size(); k++) {
                    
                    //Get vertex in A
                    int vA = verticesA_inside[k];
                    
                    //Add the displaced vertex in A
                    PointB = PointB + GNP_A.vertices[vA] + disp;
                }
                
                //Add vertices in B
                for (int k = 0; k < (int)verticesB_inside.size(); k++) {
                    
                    //Get vertex in B
                    int vB = verticesB_inside[k];
                    
                    //Add the vertex in B
                    PointB = PointB + GNP_B.vertices[vB];
                }
                
                //Get the number of vertices averaged as a double
                double n_vertices = (double)(verticesA_inside.size() + verticesB_inside.size());
                
                //Get the average
                PointB = PointB/n_vertices;
                
                //Calculate PointA by displacing B back towards A
                PointA = PointB - disp;
            }
        }
    }
    else {
        hout<<"Error in Find_point_b_for_face_in_simplex_a: simplexB has an invalid size="<<simplexB.size()<<", it size must be 4"<<endl;
        hout << "GNP for simplex of size 4 is saved in file gnp_A.vtk." << endl;
        hout << "GNP for simplex of size " << simplexB.size() << " is saved in file gnp_B.vtk." << endl;
        hout << "GNP_A=" << GNP_A.flag << " GNP_B=" << GNP_B.flag << endl;

        //Export the GNPs for debugging
        VTK_Export VTK;
        VTK.Export_single_gnp(GNP_A, "gnp_A.vtk");
        VTK.Export_single_gnp(GNP_B, "gnp_B.vtk");

        return 0;
    }
    
    
    return 1;
}
//This function finds the vertices that are contained on a face
int Hoshen_Kopelman::Get_vertices_inside_face(const vector<int> &simplexA, const vector<Edge> &edgesB, const vector<int> &normalsB, const GNP &GNP_A, const GNP &GNP_B, vector<int> &verticesA_inside)
{
    
    //Iterate over the vertices in simplexA
    for (int j = 0; j < (int)simplexA.size(); j++) {
        
        //Get current vertex of simplexA
        int vA = simplexA[j];
        
        //Variable to count the negative dot products
        int count = 0;
        
        //Iterate over the edges of simplexB and count the negative and zero dot products
        for (int i = 0; i < (int)edgesB.size(); i++) {
            
            //Get the current normal of edge i on B
            int Nb = normalsB[i];
            
            //Get vertex 1 of edge i on B
            int v1B = edgesB[i].v1;
            
            //Check if the dot product is negative or zero
            if ((GNP_A.vertices[vA] - GNP_B.vertices[v1B]).dot(GNP_B.faces[Nb].N) < Zero) {
                
                //Increase the count for vertex j in A
                count++;
            }
        }
        
        //If vertex j is inside B, add it to the vector of vertices
        if (count == 4) {
            verticesA_inside.push_back(vA);
        }
    }
    
    //Make sure at least one vertex is inside or at the boundary
    if (verticesA_inside.empty()) {
        hout<<"Error in Get_vertices_inside_face, no vertices inside face. If two faces are in contact, at least one vertex of each face should be \"inside\" the other face"<<endl;
        return 0;
    }
    
    return 1;
}
//In the case of a single vertex being inside the face of the other GNP, this function
//finds two intersections with an edge and averages them with the vertex inside the face
//Here, A is the face that contains vertex vB of B
int Hoshen_Kopelman::Average_point_with_two_interserctions(const vector<Edge> &edgesA, const vector<int> &normalsA, const GNP &gnpA, const vector<Edge> &edgesB, const vector<int> &normalsB, const GNP &gnpB, const Point_3D &disp, const int &vB, Point_3D &PB)
{
    //Variable to store the vertices of B
    int vrtB[] = {0,0};
    int it = 0;
    
    //Go through all edges in B and find the edges that intersect an edge in face A
    for (int i = 0; i < (int)edgesB.size(); i++) {
        
        //Check if edge i in B contains vertex vB
        if (edgesB[i].v1 == vB) {
            vrtB[it] = edgesB[i].v2;
            it++;
        }
        else if (edgesB[i].v2 == vB) {
            vrtB[it] = edgesB[i].v1;
            it++;
        }
    }
    
    //Variable to store the intersected edge
    int v1A = 0;
    int v2A = -1;
    
    //Iterate over the edges in A to find the intersected edge
    for (int i = 0; i < (int)edgesA.size(); i++) {
        
        //Get vertex 1 of current edge
        v1A = edgesA[i].v1;
        
        //Get the current normal
        int NA = normalsA[i];
        
        //Check the sign of the dot products
        if ((gnpB.vertices[vrtB[0]] - gnpA.vertices[v1A]).dot(gnpA.faces[NA].N) > Zero &&
            (gnpB.vertices[vrtB[1]] - gnpA.vertices[v1A]).dot(gnpA.faces[NA].N) > Zero) {
            
            //If both dot product are positive, then edge i on A
            v2A = edgesA[i].v2;
            
            //Break the loop as it is not necessary to continue searching for the edge
            break;
        }
    }
    
    //Just a check, make sure i_ed is not -1
    if (v2A == -1) {
        hout<<"Error in Average_point_with_two_interserctions: i_ed was not updated"<<endl;
        return 0;
    }
    
    //Find the two lambdas
    //Move the vertices of A towards B
    double lambda0 = Lambda_of_two_lines(PB, gnpB.vertices[vrtB[0]], gnpA.vertices[v1A] - disp, gnpA.vertices[v2A] - disp);
    double lambda1 = Lambda_of_two_lines(PB, gnpB.vertices[vrtB[1]], gnpA.vertices[v1A] - disp, gnpA.vertices[v2A] - disp);
    
    //Add the two intersection points to PB
    PB = PB + (gnpB.vertices[vrtB[0]] + (PB - gnpB.vertices[vrtB[0]])*lambda0);
    PB = PB + (gnpB.vertices[vrtB[1]] + (PB - gnpB.vertices[vrtB[1]])*lambda1);
    
    //Take the average PB
    PB = PB/3.0;
    
    return 1;
}
//Add junction to vector of junctions
int Hoshen_Kopelman::Add_gnp_junction(const vector<GNP>& gnps, const int& GNPa, const int& GNPb, const bool& p_flag, const Point_3D& PointA, const Point_3D& PointB, const double& dist, const double& tol_gnp, vector<Point_3D>& points_gnp, vector<vector<long int> >& structure_gnp)
{
    //Variables to store the GNP point numbers
    long int Pa, Pb;
    
    //Get the actual value of Pa and add PointA to vectors if not repeated
    if (!Add_gnp_point_to_vectors_if_not_repeated(PointA, tol_gnp, p_flag, Pa, points_gnp, structure_gnp[GNPa]))
    {
        hout << "Error in Add_gnp_junction when calling Add_gnp_point_to_vectors_if_not_repeated for PointA" << endl;
        return 0;
    }

    //Get the actual value of Pb and add PointB to vectors if not repeated
    if (!Add_gnp_point_to_vectors_if_not_repeated(PointB, tol_gnp, p_flag, Pb, points_gnp, structure_gnp[GNPb]))
    {
        hout << "Error in Add_gnp_junction when calling Add_gnp_point_to_vectors_if_not_repeated for PointB" << endl;
        return 0;
    }

    //Create a junction with the GNPs in contact
    Junction j(Pa, GNPa, "GNP", Pb, GNPb, "GNP", dist);

    //Add the junction to the vector of junctions
    junctions_gnp.push_back(j); 
    
    //hout << "Junction: GNPa=" << GNPa << " Pa=" << Pa << " GNPb=" << GNPb << " Pb=" << Pb << " dist=" << dist << endl;

    return 1;
}
//This function adds a GNP junction point to the vectors of GNP points and the structure vector
//However, if the point is already in the vector of points or too close to it, it is not added 
//and the existing point is used instead
int Hoshen_Kopelman::Add_gnp_point_to_vectors_if_not_repeated(const Point_3D& P, const double& tol_gnp, const bool& p_flag, long int& Pj, vector<Point_3D>& points_gnp, vector<long int>& structure_gnp)
{
    //Check if there were interpenetrations
    if (p_flag)
    {
        //Iterate over all points in GNP
        for (size_t i = 0; i < structure_gnp.size(); i++)
        {
            //Get the point number
            long int Pgnp = structure_gnp[i];

            //Check if too close to P
            //Use squared distance to reduce computational time
            if (P.squared_distance_to(points_gnp[Pgnp]) < tol_gnp)
            {
                //P and Pgnp are considered to be the same point
                //Set the point number to be Pgnp
                //This is the point number that will be used to add a junction object
                Pj = Pgnp;

                //There is no need to add a new point to the vector of points nor to the 
                //structure since the point we need already exists (as Pgnp)

                //Terminate function to avoid executing the code after the if-statement
                //that contains the for-loop
                return 1;
            }
        }
    }

    //If this part of the code is reached, then there are no interpenetrations nor were
    //repeated points found
    //Get the point number of P
    Pj = (long int)points_gnp.size();

    //Add P as a new point to the vector of points
    points_gnp.push_back(P);

    //Add Pj to the structure vector
    structure_gnp.push_back(Pj);

    return 1;
}
//This function merges CNT and GNP labels to make mixed clusters
int Hoshen_Kopelman::Make_mixed_clusters(const int &n_labels, const Cutoff_dist &cutoffs, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &sectioned_domain_cnt, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, vector<int> &labels_cnt, vector<int> &labels_gnp, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp, int &n_clusters)
{
    //Vector to store all GNP points in contact with CNTs before contact segments are compressed
    vector<Point_3D> contact_points_gnp;
    
    //Variables needed to compress contac segments
    map<long int, map<int, long int> > point_contacts;
    map<long int, map<int, double> > point_contacts_dist;
    vector<map<int, set<long int> > > gnp_cnt_point_contacts(gnps.size());
    
    //Adjacency matrix used to meger CNT and GNP clusters
    vector<set<int> > adj_labels(n_labels);
    
    //Fill the adjacency matrix and the variables needed to compress contact segments
    //hout<<"Find_adjacent_labels"<<endl;
    if (!Find_adjacent_labels(labels_cnt, labels_gnp, cutoffs, points_cnt, radii, sectioned_domain_cnt, gnps, sectioned_domain_gnp, contact_points_gnp, point_contacts, point_contacts_dist, gnp_cnt_point_contacts, adj_labels)) {
        hout<<"Error in Make_mixed_clusters when calling Find_adjacent_labels"<<endl;
        return 0;
    }
    
    //Compress contacts
    //hout<<"Compress_mixed_contacts"<<endl;
    if (!Compress_mixed_contacts(cutoffs, point_contacts, point_contacts_dist, gnp_cnt_point_contacts, contact_points_gnp, structure_gnp, points_gnp)) {
        hout<<"Error in Make_mixed_clusters when calling Compress_mixed_contacts"<<endl;
        return 0;
    }
    
    //Merge labels
    //hout<<"Merge_cnt_and_gnp_labels"<<endl;
    if (!Merge_cnt_and_gnp_labels(n_labels, labels_cnt, labels_gnp, adj_labels, n_clusters)) {
        hout<<"Error in Make_mixed_clusters when calling Merge_labels"<<endl;
        return 0;
    }
    //hout<<"Make_mixed_clusters end"<<endl;
    
    return 1;
}
//This function finds the CNTs in contact with GNPs, then an adjacency matrix is filled to
//merge the CNT and GNP clusters
int Hoshen_Kopelman::Find_adjacent_labels(
    const vector<int> &labels_cnt, 
    const vector<int> &labels_gnp, 
    const Cutoff_dist &cutoffs, 
    const vector<Point_3D> &points_cnt, 
    const vector<double> &radii, 
    const vector<vector<long int> > &sectioned_domain_cnt, 
    const vector<GNP> &gnps, 
    const vector<vector<int> > &sectioned_domain_gnp, 
    vector<Point_3D> &contact_points_gnp, 
    map<long int, map<int, long int> > &point_contacts, 
    map<long int, map<int, double> > &point_contacts_dist, 
    vector<map<int, set<long int> > > &gnp_cnt_point_contacts, 
    vector<set<int> > &adj_labels)
{
    //Create a vector of visited contacts
    vector<set<long int> > visited(gnps.size());
    
    //Use a Generate_Network object to calcualte the distance between a CNT point and a GNP
    Generate_Network GN;
    
    //Iterate over the subregions
    for (int i = 0; i < (int)sectioned_domain_gnp.size(); i++) {
        //hout<<"i="<<i<<endl;
        
        //Iterate over the GNPs in the given subregion
        for (int j = 0; j < (int)sectioned_domain_gnp[i].size(); j++) {
            //hout<<"j="<<j<<endl;
            
            //Get the current GNP
            int GNP = sectioned_domain_gnp[i][j];
            
            //Itereate over the points in the given subregion
            for (int k = 0; k < (int)sectioned_domain_cnt[i].size(); k++) {
                //hout<<"k="<<k<<endl;
                
                //Get current point
                long int P = sectioned_domain_cnt[i][k];
                
                //Get current CNT
                int CNT = points_cnt[P].flag;
                
                //Check if the GNP-CNT contact at point P has already been visited
                if (visited[GNP].find(P) != visited[GNP].end()) {
                    
                    //The GNP-CNT contact has not been visited, so calculate the distance
                    //between CNT point and GNP
                    double dist_cnt_gnp;
                    Point_3D P_gnp = GN.Get_gnp_point_closest_to_point(gnps[GNP], points_cnt[P], dist_cnt_gnp);

                    //Calculate the junction distance
                    double dist_junc = dist_cnt_gnp - radii[CNT];
                    //hout << "dist_junc=" << dist_junc <<" dist_cnt_gnp="<< dist_cnt_gnp << endl;

                    //Check that the juntion distance is at least the van der Waals distance
                    if (dist_junc < cutoffs.van_der_Waals_dist)
                    {
                        //The distance between a CNT point and a GNP might be less 
                        //than the van der Waals distance if penetrations are allowed 
                        //or Abaqus causes them
                        dist_junc = cutoffs.van_der_Waals_dist;
                    }
                    
                    //Check if there is tunneling
                    if (dist_junc - cutoffs.tunneling_dist <= Zero) {
                        
                        //Get the GNP point number
                        long int P_gnp_num = (long int)contact_points_gnp.size();
                        
                        //Save the correspoding point contacts
                        point_contacts[P][GNP] = P_gnp_num;
                        
                        //Save the junction distance
                        point_contacts_dist[P][GNP] = dist_junc;
                        
                        //Add the point to the vector of GNP points
                        contact_points_gnp.push_back(P_gnp);
                        
                        //Add the CNT point to the GNP-CNT contact
                        gnp_cnt_point_contacts[GNP][CNT].insert(P);
                        
                        //Get the CNT and GNP labels
                        int l_cnt = labels_cnt[CNT];
                        int l_gnp = labels_gnp[GNP];
                        
                        //Fill the adjacency matrix for the labels
                        adj_labels[l_cnt].insert(l_gnp);
                        adj_labels[l_gnp].insert(l_cnt);
                        
                        //Add the contact to the vector of visited contacts
                        visited[GNP].insert(P);
                    }
                }
            }
        }
    }
    
    return 1;
}
//This function compresses the "contact segments" for CNT-CNT contacts
int Hoshen_Kopelman::Compress_mixed_contacts(const Cutoff_dist &cutoffs, map<long int, map<int, long int> > &point_contacts, map<long int, map<int, double> > &point_contacts_dist, vector<map<int, set<long int> > > &gnp_cnt_point_contacts, const vector<Point_3D> &contact_points_gnp, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp)
{
    //Iterate over all GNPs
    for (int i = 0; i < (int)gnp_cnt_point_contacts.size(); i++) {
        
        //Just for clarity in the following, create a variable for the current GNP
        int GNPi = i;
        
        //Check if GNPi has contacts with any CNT
        if (!gnp_cnt_point_contacts[GNPi].empty()) {
            
            //Iterate over all CNTs in contact with GNPi
            for (map<int, set<long int> >::const_iterator i_map = gnp_cnt_point_contacts[GNPi].begin(); i_map != gnp_cnt_point_contacts[GNPi].end(); i_map++) {
                
                //Get the CNT number in contact with GNPi
                int CNTj = i_map->first;
                
                //Start the iterator at the beginning of the set
                set<long int>::const_iterator i_set = i_map->second.begin();
                
                //Get the initial point on the CNT in contact with the GNP
                long int Pj0 = *i_set;
                
                //Get the intial point on the GNP in contact with the CNT
                long int Pi0 = point_contacts.at(Pj0).at(GNPi);
                
                //Get the initial junction distance
                long int d_junction = point_contacts_dist.at(Pj0).at(GNPi);
                
                //Variables to store the "compress" contact segment
                long int Pi_junc = Pi0, Pj_junc = Pj0;
                double d_junc_min = d_junction;
                
                //Iterate over the points in CNTj in contact with GNPi and compress segements
                for (i_set++; i_set != i_map->second.end(); i_set++) {
                    
                    //Get the point in CNTj
                    long int Pj1 = *i_set;
                    
                    //Get the point in GNPi
                    long int Pi1 = point_contacts.at(Pj1).at(GNPi);
                    
                    //Get the junction distance
                    d_junction = point_contacts_dist.at(Pj1).at(GNPi);
                    
                    //Check if the points are close enough to be in the same segment
                    if (Pj1 - Pj0 >= cutoffs.min_points) {
                        
                        //A new segment needs to be started
                        //But first, "compress" the previous segment
                        
                        //Get the actual GNP point number
                        long int P_gnp_num = (long int)points_gnp.size();
                        
                        //Add the GNP point from the vector of GNP contact points
                        points_gnp.push_back(contact_points_gnp[Pi_junc]);
                        
                        //Set the flag of the point
                        points_gnp.back().flag = GNPi;
                        
                        //Add the GNP point number to the structure
                        structure_gnp[GNPi].push_back(P_gnp_num);
                        
                        //Add the shortest junction found to the vector of junctions
                        Junction j(Pj_junc, CNTj, "CNT", P_gnp_num, GNPi, "GNP", d_junc_min);
                        junctions_mixed.push_back(j);
                        
                        //Add the junction point on the CNT to the vectors of elements
                        elements_cnt[CNTj].insert(Pj_junc);
                        
                        //Start the new segment by resetting the variables to the values
                        //of the junction that starts the new segment
                        d_junc_min = d_junction;
                        Pi_junc = Pi1;
                        Pj_junc = Pj1;
                    }
                    else {
                        
                        //Current contact is part of the same segment
                        //Check if junction distance and points in contact need updating
                        //Keep the points with the shortest junction distance since this
                        //is the path of lowest electrical resistance
                        if (d_junction < d_junc_min) {
                            
                            //A shorter distance was found, so update the variables for
                            //the junction for the compressed contact segment
                            d_junc_min = d_junction;
                            Pi_junc = Pi1;
                            Pj_junc = Pj1;
                        }
                    }
                }
                
                //"Compress" the last (or only) segment
                
                //Get the actual GNP point number
                long int P_gnp_num = (long int)points_gnp.size();
                
                //Add the GNP point from the vector of GNP contact points
                points_gnp.push_back(contact_points_gnp[Pi_junc]);
                
                //Set the flag of the point
                points_gnp.back().flag = GNPi;
                
                //Add the shortest junction found in the last (or only) segment
                //to the vector of junctions
                Junction j(Pj_junc, CNTj, "CNT", P_gnp_num, GNPi, "GNP", d_junc_min);
                junctions_mixed.push_back(j);
                
                //Add the junction point on the CNT to the vectors of elements
                elements_cnt[CNTj].insert(Pj_junc);
            }
        }
    }
    
    return 1;
}
//This function merges the labels using DFS and then renumbers the CNT and GNP labels with
//the merged label numbers
int Hoshen_Kopelman::Merge_cnt_and_gnp_labels(const int &n_labels, vector<int> &labels_cnt, vector<int> &labels_gnp, const vector<set<int> > &adj_labels, int &n_clusters)
{
    //Label map for merging CNT and GNP labels
    vector<int> label_map(n_labels);
    
    //Perform a DFS on the adjacency matrix of the labels
    vector<vector<int> > mixed_clusters;
    //hout<<"DFS_on_labels"<<endl;
    if (!DFS_on_labels(n_labels, adj_labels, mixed_clusters)) {
        hout<<"Error in Merge_labels when calling DFS_on_labels"<<endl;
        return 0;
    }
    
    //The number of clusters is equal to the size of mixed_clusters
    n_clusters = (int)mixed_clusters.size();
    
    //Fill the map for the merged labels
    for (int i = 0; i < (int)mixed_clusters.size(); i++) {
        
        //Iterate over all labels in mixed cluster i
        for (int j = 0; j < (int)mixed_clusters[i].size(); j++) {
            
            //Get the current label of cluster i
            int label = mixed_clusters[i][j];
            
            //Set the label map as label->i
            label_map[label] = i;
        }
    }
    
    //Renumber the CNT labels with the merged labels using label_map
    //hout<<"Map_labels CNT"<<endl;
    if (!Map_labels(label_map, labels_cnt)) {
        hout<<"Error in Merge_labels when calling Map_labels (CNTs)"<<endl;
        return 0;
    }
    
    //Renumber the GNP labels with the merged labels using label_map
    //hout<<"Map_labels GNP"<<endl;
    if (!Map_labels(label_map, labels_gnp)) {
        hout<<"Error in Merge_labels when calling Map_labels (GNPs)"<<endl;
        return 0;
    }
    
    return 1;
}
//This function performs DFS on the adjacency matrix of labels to group the CNT and GNP labels
//that are part of the same cluster of mixed particles
int Hoshen_Kopelman::DFS_on_labels(const int &n_labels, const vector<set<int> > &adj_labels, vector<vector<int> > &mixed_clusters)
{
    //Set up the vector of visited labels fo the DFS
    vector<int> visited_labels(n_labels, 0);
    
    //Iterate over the labels
    for (int i = 0; i < n_labels; i++) {
        
        //Check if label i has already been visited
        if (!visited_labels[i]) {
            
            //Set label i as visited
            visited_labels[i] = 1;
            
            //Increase the size of the vector of mixed clusters to start a new cluster
            mixed_clusters.push_back(vector<int>());
            
            //Add label i to the new cluster
            mixed_clusters.back().push_back(i);
            
            //Explore all labels connected to label i
            if (!Explore_labels(i, adj_labels, mixed_clusters, visited_labels)) {
                hout<<"Error in Explore_labels when calling Explore_labels recursively"<<endl;
                return 0;
            }
        }
    }
    
    return 1;
}
int Hoshen_Kopelman::Explore_labels(const int &label, const vector<set<int> > &adj_labels, vector<vector<int> > &mixed_clusters, vector<int> &visited_labels)
{
    //Go through the labels connected to "label"
    for (set<int>::const_iterator it = adj_labels[label].begin(); it != adj_labels[label].end(); it++) {
        
        //Get current label
        int new_label = *it;
        
        //Check if the connected label has been visited
        if (!visited_labels[new_label]) {
            
            //The label has not been visited, so set it as visited
            visited_labels[new_label] = 1;
            
            //new_label is part of the current mixed cluster, so add it to the last mixed cluster
            mixed_clusters.back().push_back(new_label);
            
            //Explore the labels connected to new_label
            if (!Explore_labels(new_label, adj_labels, mixed_clusters, visited_labels)) {
                hout<<"Error in Explore_labels when calling Explore_labels recursively"<<endl;
                return 0;
            }
        }
    }
    
    return 1;
}
int Hoshen_Kopelman::Map_labels(const vector<int> label_map, vector<int> &labels)
{
    //Renumber the labels with the merged labels using label_map
    //hout<<"label_map.size="<<label_map.size()<<" labels.size="<<labels.size()<<endl;
    for (int i = 0; i < (int)labels.size(); i++) {
        
        //Check if CNT i has a label assigned
        if (labels[i] != -1) {
            
            //For clarity, get the old label of CNT i
            int old_label = labels[i];
            //hout<<"old_label="<<old_label<<endl;
            
            //Get the new label number using label_map
            int new_label = label_map[old_label];
            
            //Assign the new label to the CNT
            labels[i] = new_label;
        }
    }
    
    return 1;
}
//In this function the clusters are made using the labels assugned by HK76
int Hoshen_Kopelman::Make_particle_clusters(const int &n_clusters, const vector<int> &particles_inside, vector<int> &labels, vector<vector<int> > &isolated, vector<vector<int> > &clusters_particles)
{
    //Assing the correct size to the vector of clusters and the vector of isolated particles:
    //the number of clusters is the same as the number of proper labels
    clusters_particles.assign(n_clusters, vector<int>());
    isolated.push_back(vector<int>());
    
    //Now scan the vector of particles_inside. Check the label of each particle to make the clusters.
    for (int i = 0; i < (int)particles_inside.size(); i++) {
        
        //Get the current particle
        int particle = particles_inside[i];
        
        //Get the label of the current particle
        int L = labels[particle];
        
        //If a label L is -1, it means that particle is isolated
        //At this point "isolated" means that the particle is not part of any cluster
        //Only when L is diffenrent form -1, particle_i belongs to a cluster
        if (L != -1){
            
            //At this point, labels are consecutive and starting from 0
            //and thus the label indicates the cluster number
            clusters_particles[L].push_back(particle);
        } else {
            
            //Add isolated CNTs to the last elment in the vector of isolated particles
            isolated.back().push_back(particle);
        }
    }
    
    return 1;
}
//This function determines percolation and, in case there is percolation,
//assigns the family number to each cluster
int Hoshen_Kopelman::Find_percolated_clusters(const int &n_clusters, const vector<vector<int> > &boundary_cnt, const vector<vector<int> > &boundary_gnp, const vector<int> &labels_cnt, const vector<int> &labels_gnp, const vector<vector<int> > &clusters_cnt_tmp, const vector<vector<int> > &clusters_gnp_tmp, map<int,int> &percolated_labels)
{
    //Array of opposite boundaries
    //{2,4}: x-boundaries
    //{3,5}: y-boundaries
    //{0,1}: z-boundaries
    int boundary_pairs[][2] = {{2,4},{3,5},{0,1}};
    
    //Vector of percolated directions
    vector<set<int> > percolated_dirs(n_clusters, set<int>());
    
    //Find the clusters that are connected to boundaries and if they connect
    //opposite boundaries, i.e., clusters that percolate
    //hout<<"Find_clusters_connected_to_boundaries"<<endl;
    if (!Find_clusters_connected_to_boundaries(boundary_cnt, boundary_gnp, labels_cnt, labels_gnp, boundary_pairs, percolated_dirs)) {
        hout<<"Error in Find_percolated_clusters when calling Find_clusters_connected_to_boundaries"<<endl;
        return 0;
    }
    
    //Add the percolated clusters to the vectors of percolated clusters, and the non-percoalted
    //clusters to the vectors of isolated particles
    //Also, determine the family of each percoalted cluster
    //hout<<"Determine_family_of_percolated_clusters"<<endl;
    if (!Determine_family_of_percolated_clusters(n_clusters, clusters_cnt_tmp, clusters_gnp_tmp, percolated_dirs, percolated_labels)) {
        hout<<"Error in Find_percolated_clusters when calling Determine_family_of_percolated_clusters"<<endl;
        return 0;
    }
    
    return 1;
}
//This function finds the clusters that are connected to a boundary
//It is also determined if a cluster connects two opposite boundaries
int Hoshen_Kopelman::Find_clusters_connected_to_boundaries(const vector<vector<int> > &boundary_cnt, const vector<vector<int> > &boundary_gnp, const vector<int> &labels_cnt, const vector<int> &labels_gnp, const int boundary_pairs[][2], vector<set<int> > &percolated_dirs)
{
    //Iterate over the three boundary pairs
    for (int i = 0; i < 3; i++) {
        
        //Get one boundary
        int b1 = boundary_pairs[i][0];
        
        //Set to store the cluster numbers (i.e., labels) that are connected
        //to boundary b1
        set<int> boundary1;
        
        //Check if there are CNTs
        if (boundary_cnt.size()) {
            
            //Add the CNT clusters connected to boundary b1 to the set boundary1
            //hout<<"Add_clusters_in_boundary dir="<<i<<" CNT boundary_cnt[b1].size="<<boundary_cnt[b1].size()<<endl;
            if (!Add_clusters_in_boundary(boundary_cnt[b1], labels_cnt, boundary1)) {
                hout<<"Error in Find_clusters_connected_to_boundaries when calling Add_clusters_in_boundary (CNT)"<<endl;
                return 0;
            }
        }
        
        //Check if there are GNPs
        if (boundary_gnp.size()) {
            
            //Add the GNP clusters connected to boundary b1 to the set boundary1
            //hout<<"Add_clusters_in_boundary dir="<<i<<" GNP boundary_gnp[b1="<<b1<<"].size="<<boundary_gnp[b1].size()<<endl;
            if (!Add_clusters_in_boundary(boundary_gnp[b1], labels_gnp, boundary1)) {
                hout<<"Error in Find_clusters_connected_to_boundaries when calling Add_clusters_in_boundary (GNP)"<<endl;
                return 0;
            }
        }
        
        //Check if there were any clusters connected to boundary b1
        if (!boundary1.empty()) {
            
            //There are clusters connected to boundary1, so there might be percolation
            
            //Get the boundary opposite to b1
            int b2 = boundary_pairs[i][1];
            
            //Check if there are CNTs
            if (boundary_cnt.size()) {
                
                //Find the CNT clusters connected to boundary b2 and add a percolated direction
                //if percolation is determined
                //hout<<"Add_percolated_direction dir="<<i<<" CNT boundary_cnt[b2].size="<<boundary_cnt[b2].size()<<endl;
                if (!Add_percolated_direction(i, boundary_cnt[b2], labels_cnt, boundary1, percolated_dirs)) {
                    hout<<"Error in Find_clusters_connected_to_boundaries when calling Add_percolated_direction (CNT)"<<endl;
                    return 0;
                }
            }
            
            //Check if there are GNPs
            if (boundary_gnp.size()) {
                
                //Find the GNP clusters connected to boundary b2 and add a percolated direction
                //if percolation is determined
                //hout<<"Add_percolated_direction dir="<<i<<" GNP boundary_gnp[b2="<<b2<<"].size="<<boundary_gnp[b2].size()<<endl;
                if (!Add_percolated_direction(i, boundary_gnp[b2], labels_gnp, boundary1, percolated_dirs)) {
                    hout<<"Error in Find_clusters_connected_to_boundaries when calling Add_percolated_direction (GNP)"<<endl;
                    return 0;
                }
            }
        }
    }
    
    return 1;
}
//This function adds the clusters connected to a boundary (as given by the vector
//boundary_particles) into the set boundary1
int Hoshen_Kopelman::Add_clusters_in_boundary(const vector<int> &boundary_particles, const vector<int> &labels, set<int> &boundary1)
{
    //Iterate over the particles on the boundary vector (if any)
    for (int j = 0; j < (int)boundary_particles.size(); j++) {
        
        //Get current particle on boundary vector
        int part = boundary_particles[j];
        
        //Get the cluster number, i.e., label, of particle
        int L = labels[part];
        //hout<<"Particle="<<part<<" L="<<L<<endl;
        
        //If the particle has a label assigned, then add it to the set of clusters
        //connected to boundary b1
        if (L != -1) {
            boundary1.insert(L);
        }
    }
    
    return 1;
}
//This function finds the clusters connected to a boundary as given by the boundary
//vector "boundary_particles"
//These clusters are compared with the clusters connected to the opposite boundary,
//as indicated by the set boundary1
//When a cluster is connected to both boundaries, then it percolated along d and this is
//added to the vector of percoalted directions (percolated_dirs)
int Hoshen_Kopelman::Add_percolated_direction(const int &d, const vector<int> &boundary_particles, const vector<int> &labels, const set<int> &boundary1, vector<set<int> > &percolated_dirs)
{
    //Iterate over the particles on the boundary vector (if any)
    for (int j = 0; j < (int)boundary_particles.size(); j++) {
        
        //Get current particle on boundary b
        int part = boundary_particles[j];
        
        //Get the cluster number, i.e., label, of particle
        int L = labels[part];
        //hout<<"Particle="<<part<<" L="<<L<<endl;
        
        //Check if cluster L is also in boundary b1
        if (boundary1.find(L) != boundary1.end()) {
            
            //Cluster L percolates along direction d
            
            //Map direction d
            int map = 10*d + 1;
            
            //Add the mapped value of the direction to the vector of percolated direction
            percolated_dirs[L].insert(map);
        }
    }
    return 1;
}
//This function determines the family of a percoalted cluster
int Hoshen_Kopelman::Determine_family_of_percolated_clusters(const int &n_clusters, const vector<vector<int> > &clusters_cnt_tmp, const vector<vector<int> > &clusters_gnp_tmp, const vector<set<int> > &percolated_dirs, map<int,int> &percolated_labels)
{
    //Variable to count the percolated clusters
    int percolated = 0;
    
    //Go though the vector of percolated directions and determine the family of each cluster
    for (int i = 0; i < n_clusters; i++) {
        
        //Check if cluster i percolates
        //hout<<"percolated_dirs[i].size="<<percolated_dirs[i].size()<<endl;
        if (percolated_dirs[i].empty()) {
            
            //Cluster i does not percolate, so add it to the vector of isolated particles
            if (clusters_cnt_tmp.size()) {
                isolated_cnt.push_back(clusters_cnt_tmp[i]);
            }
            if (clusters_gnp_tmp.size()) {
                isolated_gnp.push_back(clusters_gnp_tmp[i]);
            }
        }
        else {
            
            //Cluster i percolates, so add a map of label i to percolated
            percolated_labels[i] = percolated;
            
            //Calculate the sum of mapped directions
            int sum = 0;
            for (set<int>::iterator it = percolated_dirs[i].begin(); it != percolated_dirs[i].end(); it++) {
                
                //Add the set element to the sum
                sum = sum + *it;
            }
            
            //Variable to store the family number
            int fam = -1;
            
            //Get the family number from the sum
            //hout<<"Family_map"<<endl;
            if (!Family_map(sum, fam)) {
                hout<<"Error in Find_percolated_clusters when calling Family_map"<<endl;
                return 0;
            }
            
            //Check a valid family was found
            if (fam == -1) {
                hout<<"Error in Find_percolated_clusters. No valid family was found. "<<endl;
                return 0;
            }
            //Add the family to the vector of families
            family.push_back(fam);
            
            //Add the cluster to the vectors of percolated clusters
            if (clusters_cnt_tmp.size()) {
                clusters_cnt.push_back(clusters_cnt_tmp[i]);
            }
            if (clusters_gnp_tmp.size()) {
                clusters_gnp.push_back(clusters_gnp_tmp[i]);
            }
            
            //Increase the count of percolated clusters
            percolated++;
        }
    }
    
    return 1;
}
//This function maps the sum of directions to the integer that represents the percolated family
int Hoshen_Kopelman::Family_map(const int &perc_sum, int &family)
{
    //Return the family number according to the sum
    switch (perc_sum) {
        case 1:
            family = 0;
            break;
        case 11:
            family = 1;
            break;
        case 21:
            family = 2;
            break;
        case 12:
            family = 3;
            break;
        case 22:
            family = 4;
            break;
        case 32:
            family = 5;
            break;
        case 33:
            family = 6;
            break;
            
        default:
            hout<<"Unknown sum of percolated directions: "<<perc_sum<<endl;
            return 0;
    }
    return 1;
}
//This function groups the junctions into the clusters they belong to
int Hoshen_Kopelman::Group_junctions(const vector<Point_3D> &points_cnt, const vector<Point_3D> &points_gnp, const vector<int> &labels_cnt, const vector<int> &labels_gnp, const map<int,int> &percolated_labels)
{
    //CNT-CNT junction
    if (junctions_cnt.size()) {
        
        //Set the vector of clusters for junctions to the correct size
        cluster_cnt_junctions.assign(clusters_cnt.size(), vector<int>());
        
        //Add junctions to cluster
        //hout<<"Group_junctions_same_particle CNT"<<endl;
        if (!Group_junctions_same_particle(points_cnt, labels_cnt, junctions_cnt, percolated_labels, cluster_cnt_junctions)) {
            hout<<"Error in Group_junctions when calling Group_junctions_same_particle (CNTs)"<<endl;
            return 0;
        }
    }
    
    //GNP-GNP junctions
    if (junctions_gnp.size()) {
        
        //Set the vector of clusters for junctions to the correct size
        cluster_gnp_junctions.assign(clusters_gnp.size(), vector<int>());
        
        //Add junctions to cluster
        //hout<<"Group_junctions_same_particle GNP"<<endl;
        if (!Group_junctions_same_particle(points_gnp, labels_gnp, junctions_gnp, percolated_labels, cluster_gnp_junctions)) {
            hout<<"Error in Group_junctions when calling Group_junctions_same_particle (GNPs)"<<endl;
            return 0;
        }
    }
    
    //CNT-GNP junctions
    if (junctions_mixed.size()) {
        
        //Set the vector of clusters for junctions to the correct size
        cluster_mix_junctions.assign(clusters_cnt.size(), vector<int>());
        
        //Add junctions to cluster
        //hout<<"Group_junctions_mix_particle"<<endl;
        if (!Group_junctions_mix_particle(points_cnt, points_gnp, labels_cnt, labels_gnp, percolated_labels)) {
            hout<<"Error in Group_junctions when calling Group_junctions_mix_particle"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function groups the junctions into the clusters they belong to when the junction is
//between two particles of the same type, i.e., CNT-CNT and GNP-GNP junctions
int Hoshen_Kopelman::Group_junctions_same_particle(const vector<Point_3D> &points, const vector<int> &labels, const vector<Junction> &junctions, const map<int,int> &percolated_labels, vector<vector<int> > &cluster_junction)
{
    //Iterate over all junctions
    for (int i = 0; i < (int)junctions.size(); i++) {
        
        //Get the Particle number of the two particles
        int Pa1 = junctions[i].N1;
        int Pa2 = junctions[i].N2;
        
        //Get the cluster numbers of the two particles
        int L1 = labels[Pa1];
        int L2 = labels[Pa2];
        //hout<<"L1="<<L1<<" L2="<<L2<<endl;
        
        //Check that the two particles belong to the same cluster
        if (L1 != L2) {
            hout<<"Error in Group_junctions_same_particle: a junction between particles in different clusters was found. This cannot happen as two particles that have a junction belong to the same cluster"<<endl;
            hout<<"Particle 1:"<<junctions[i].type1<<" number:"<<Pa1<<" point:"<<junctions[i].P1<<" coordinates:"<<points[junctions[i].P1].str()<<endl;
            hout<<"Particle 2:"<<junctions[i].type2<<" number:"<<Pa2<<" point:"<<junctions[i].P2<<" coordinates:"<<points[junctions[i].P2].str()<<endl;
            return 0;
        }
        
        //Check if the junction is at a percolated clusters
        if (percolated_labels.find(L1) != percolated_labels.end()) {
            
            //Add junction number to corresponding cluster
            //hout<<"junction i="<<i<<" with label L="<<L1<<" added with percolated label="<<percolated_labels.at(L1)<<endl;
            cluster_junction[percolated_labels.at(L1)].push_back(i);
        }
    }
    
    return 1;
}
//This function groups the junctions into the clusters they belong to when the junction is
//between two particles of different type, i.e., CNT-GNP junctions
int Hoshen_Kopelman::Group_junctions_mix_particle(const vector<Point_3D> &points_cnt, const vector<Point_3D> &points_gnp, const vector<int> &labels_cnt, const vector<int> &labels_gnp, const map<int,int> &percolated_labels)
{
    //Iterate over all mixed junctions
    for (int i = 0; i < (int)junctions_mixed.size(); i++) {
        
        //Get the Particle number of the two particles
        int Pa1 = junctions_mixed[i].N1;
        int Pa2 = junctions_mixed[i].N2;
        
        //Get the cluster numbers of the two particles
        int L1 = labels_cnt[Pa1];
        int L2 = labels_gnp[Pa2];
        
        //Check that the two particles belong to the same cluster
        if (L1 != L2) {
            hout<<"Error in Group_junctions_mix_particle: a junction between particles in different clusters was found. This cannot happen as two particles that have a junction belong to the same cluster"<<endl;
            hout<<"CNT:"<<Pa1<<" point:"<<junctions_mixed[i].P1<<" coordinates:"<<points_cnt[junctions_mixed[i].P1].str()<<endl;
            hout<<"GNP:"<<Pa2<<" number:"<<Pa2<<" point:"<<junctions_mixed[i].P2<<" coordinates:"<<points_gnp[junctions_mixed[i].P2].str()<<endl;
            return 0;
        }
        
        //Check if the junction is at a percolated clusters
        if (percolated_labels.find(L1) != percolated_labels.end()) {
            
            //Add junction number to corresponding cluster
            cluster_mix_junctions[percolated_labels.at(L1)].push_back(i);
        }
    }
    
    return 1;
}
//-------------------------------------------------------
//HK76
//Function for the Hoshen-Kopelman (HK76) algorithm only
//It is assumed that this function is used only when a contact is found. Otherwise the results will be wrong and probably one cluster with all CNTs will be generated
int Hoshen_Kopelman::HK76(const int &CNT1, const int &CNT2, int &new_label, vector<int> &labels,  vector<int> &labels_labels) {
    //L is just a label variable and new_label will take the value of the newest cluster
    int L;
    
    //If both labels of a CNT are -1, then both CNT will form a new cluster labeled with the current value in new_label
    //If only one of the CNTs has label -1, then use the same label as the other CNT
    //If CNTs have different label and both are not -1, then the labels need to merge
    //labels_labels[i] keeps the size of the label i. If two labels are merged, then labels_labels[i] will have a negative vaule equal to the negative of the smallest merged label
    if ( (labels[CNT1] == -1) && (labels[CNT2] == -1) ) {
        labels[CNT1] = new_label;
        labels[CNT2] = new_label;
        new_label++;
        labels_labels.push_back(2);
    } else if (labels[CNT1] == -1) {
        //the proper label is in labels[CNT2]
        L = labels[CNT2];
        if (labels_labels[L] <= 0) {
            //hout << "L=" << L << " LL[L]=" << labels_labels[L] << endl;
            L = Find_root(L, labels_labels);
            //hout << "L_root=" << L << " LL[L_root]=" << labels_labels[L] << endl;
        }
        labels[CNT1] = L;
        labels_labels[L] = labels_labels[L] + 1;
    } else if (labels[CNT2] == -1) {
        //the proper label is in labels[CNT1]
        L = labels[CNT1];
        if (labels_labels[L] <= 0) {
            //hout << "L=" << L << " LL[L]=" << labels_labels[L] << endl;
            L = Find_root(L, labels_labels);
            //hout << "L_root=" << L << " LL[L_root]=" << labels_labels[L] << endl;
        }
        labels[CNT2] = L;
        labels_labels[L]= labels_labels[L] + 1;
    } else if (labels[CNT1] != labels[CNT2]) {
        //Solve label "conflict" and merge labels
        int L1 = Find_root(labels[CNT1], labels_labels);
        int L2 = Find_root(labels[CNT2], labels_labels);
        if (!Merge_labels(L1, L2, labels_labels)){
            hout << "Error in HK76" << endl;
            return 0;
        }
    }
    return 1;
}

//Find in the label L is a root or it points to another label
//If it points to another label, find that label (i.e. the root)
int Hoshen_Kopelman::Find_root(const int &L, vector<int> &labels_labels)
{
    int L_root = L;
    while (labels_labels[L_root] <= 0){
        L_root = -labels_labels[L_root];
        //If labels_labels[L_root] = 0, then the root is zero, not necesarily L
        if (labels_labels[L_root] == 0)
            return 0;
    }
    
    return L_root;
}

//this function merges two clusters
int Hoshen_Kopelman::Merge_labels(const int &root1, const int &root2, vector<int> &labels_labels)
{
    if ( (labels_labels[root1] <= 0) || (labels_labels[root2] <= 0) ) {
        hout << "Error on merging clusters. Both labels are negative, at this point this should not happen."<<endl;
        hout << "root1=" << root1 << " LL[root1]=" << labels_labels[root1];
        hout << " root2=" << root2 << " LL[root2]=" << labels_labels[root2]<< endl;
        return 0;
    }
    if (root1 < root2) {
        //In this case the root is root1
        labels_labels[root1] = labels_labels[root1] + labels_labels[root2];
        labels_labels[root2] = -root1;
    } else if (root2 < root1) {
        //In this case the root is root2;
        labels_labels[root2] = labels_labels[root2] + labels_labels[root1];
        labels_labels[root1] = -root2;
    }
    //If root1 is equal to root2, the CNTs are already in the same cluster so there is nothing to do
    
    return 1;
}
//-------------------------------------------------------
//Visualization files
//
//This function is used to export clusters before determining percolation
int Hoshen_Kopelman::Export_clusters(const int &percolation, const int &iter, const vector<vector<int> > &clusters_cnt, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<vector<int> > &clusters_gnp, const vector<GNP> &gnps)
{
    //VTK export object
    VTK_Export VTK_E;
    
    //Set the prefix of the filename for clusters
    string pref_cl = (percolation)? "perc_cluster_": "cluster_";
    
    //Set the prefix of filename for isolated particles
    string pref_iso = (percolation)? "isolated_and_non_percolated_": "isolated_";
    
    //Export CNT clusters
    for (int i = 0; i < (int)clusters_cnt.size(); i++) {
        
        //Prepare filename for cluster i
        string filename = pref_cl + "cnt_" + to_string(i) + ".vtk";
        
        //Export cluster i with speficied name
        if (!VTK_E.Export_cnts_in_cluster(points_cnt, structure_cnt, clusters_cnt[i], filename)) {
            hout<<"Error in Export_clusters when calling VTK_E.Export_cnts_in_cluster i="<<i<<endl;
            return 0;
        }
    }
    
    //Get all isolated CNTs into a single cluster
    vector<int> isolated_cnts_all;
    if (!Combine_into_one_cluster(isolated_cnt, isolated_cnts_all)) {
        hout<<"Error in Export_clusters when calling Combine_into_one_cluster (CNTs)"<<endl;
        return 0;
    }
    
    //Get name for isolated CNTs
    string filename_iso = pref_iso + "cnt.vtk";
    
    //Export isolated CNTs, if any
    if (isolated_cnts_all.size()) {
        if (!VTK_E.Export_cnts_in_cluster(points_cnt, structure_cnt, isolated_cnts_all, filename_iso)) {
            hout<<"Error in Export_clusters when calling VTK_E.Export_cnts_in_cluster (isolated)"<<endl;
            return 0;
        }
    }
    
    //Export GNP clusters
    for (int i = 0; i < (int)clusters_gnp.size(); i++) {
        
        //Prepare filename for cluster i
        string filename = pref_cl + "gnp_" + to_string(i) + ".vtk";
        
        //Export cluster i with speficied name
        if (!VTK_E.Export_gnps_in_cluster(gnps, clusters_gnp[i], filename)) {
            hout<<"Error in Export_clusters when calling VTK_E.Export_gnps_in_cluster i="<<i<<endl;
            return 0;
        }
    }
    
    //Get all isolated GNPs into a single cluster
    vector<int> isolated_gnps_all;
    if (!Combine_into_one_cluster(isolated_gnp, isolated_gnps_all)) {
        hout<<"Error in Export_clusters when calling Combine_into_one_cluster (GNPs)"<<endl;
        return 0;
    }
    
    //Get name for isolated CNTs
    filename_iso = pref_iso + "gnp.vtk";
    
    //Export isolated GNPs if any
    if (isolated_gnps_all.size()) {
        if (!VTK_E.Export_gnps_in_cluster(gnps, isolated_gnps_all, filename_iso)) {
            hout<<"Error in Export_clusters when calling VTK_E.Export_gnps_in_cluster (isolated)"<<endl;
            return 0;
        }
    }
    
    //Move all visualization files to the folder of the iteration
    char command[100];
    if (percolation) {
        sprintf(command, "mkdir percolated_%.4d", iter);
        system(command);
        sprintf(command, "mv isolated_*.vtk percolated_%.4d", iter);
        system(command);
        
        //Check if there are percolated clusters and, thus, the corresponding visualization files
        if (clusters_cnt.size() || clusters_gnp.size()) {
            sprintf(command, "mv perc_cluster_*.vtk percolated_%.4d", iter);
            system(command);
        }
    }
    else {
        sprintf(command, "mkdir clusters_%.4d", iter);
        system(command);
        sprintf(command, "mv isolated_*.vtk clusters_%.4d", iter);
        system(command);
        
        //Check if there are clusters and, thus, the corresponding visualization files
        if (clusters_cnt.size() || clusters_gnp.size()) {
            sprintf(command, "mv cluster_*.vtk clusters_%.4d", iter);
            system(command);
        }
    }
    
    return 1;
}
int Hoshen_Kopelman::Combine_into_one_cluster(const vector<vector<int> > &clusters, vector<int> &cluster)
{
    //Iterate over all vectors in clusters
    for (int i = 0; i < (int)clusters.size(); i++) {
        
        //Add all elements in clusters[i] into cluster
        cluster.insert(cluster.end(), clusters[i].begin(), clusters[i].end());
    }
    
    return 1;
}
